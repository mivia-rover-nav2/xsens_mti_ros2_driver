#include "host_clock_sync.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace
{
constexpr double kTickS = 1e-4;              // SampleTimeFine runs at 10 kHz
constexpr int64_t kMaxBackwardTicks = 1000;  // 0.1 s: more than reordering, so the MTi restarted
}

HostClockSync::HostClockSync() : HostClockSync(Config()) {}

HostClockSync::HostClockSync(const Config &config) : m_config(config), m_resets(0)
{
	if (m_config.window_s <= 0.0)
		m_config.window_s = 1.0;
	if (m_config.history < 2)
		m_config.history = 2;
	clear();
}

void HostClockSync::clear()
{
	m_started = false;
	m_prevTicks = 0;
	m_ticks = 0;
	m_hostRefNs = 0;
	m_window = 0;
	m_windowHasPoint = false;
	m_windowMin = Point{0.0, 0.0};
	m_aboveLine = false;
	m_aboveSince = 0.0;
	m_delaySum = 0.0;
	m_delayCount = 0;
	m_meanDelay = 0.0;
	m_minima.clear();
	m_a = 0.0;
	m_b = 0.0;
	m_slewStart = 0.0;
	m_slewAmount = 0.0;
	m_windows = 0;
}

void HostClockSync::reset()
{
	clear();
	++m_resets;
}

int64_t HostClockSync::update(uint32_t sample_time_fine, int64_t host_ns)
{
	if (m_started)
	{
		// Unsigned subtraction wraps modulo 2^32, so the cast gives the signed
		// step across a counter rollover as well.
		const int64_t step = static_cast<int32_t>(sample_time_fine - m_prevTicks);
		if (step < -kMaxBackwardTicks)
			reset();  // the MTi restarted its counter
		else
		{
			m_ticks += step;
			m_prevTicks = sample_time_fine;
		}
	}
	if (!m_started)
		start(sample_time_fine, host_ns);

	double sensor_s = static_cast<double>(m_ticks) * kTickS;
	double offset_s = static_cast<double>(host_ns - m_hostRefNs) * 1e-9 - sensor_s;

	if (m_windows > 0)
	{
		const double residual = offset_s - predictOffset(sensor_s);
		bool stepped = false;
		if (residual < -m_config.step_threshold_s)
		{
			// Arrived well before the fitted line: impossible unless the host
			// clock stepped back (or the sensor's forward).
			stepped = true;
		}
		else if (residual > m_config.step_threshold_s)
		{
			// Late packets alone (a USB stall) recover within a burst; staying
			// far above the line means the host clock stepped forward (or the
			// sensor's back).
			if (!m_aboveLine)
			{
				m_aboveLine = true;
				m_aboveSince = sensor_s;
			}
			stepped = sensor_s - m_aboveSince >= m_config.step_confirm_s;
		}
		else
		{
			m_aboveLine = false;
			m_delaySum += residual;
			++m_delayCount;
		}
		if (stepped)
		{
			reset();
			start(sample_time_fine, host_ns);
			sensor_s = 0.0;
			offset_s = 0.0;
		}
	}

	const long window = static_cast<long>(std::floor(sensor_s / m_config.window_s));
	if (m_windowHasPoint && window > m_window)
		closeWindow(sensor_s);
	if (!m_windowHasPoint)
	{
		m_window = static_cast<long>(std::floor(sensor_s / m_config.window_s));
		m_windowMin = Point{sensor_s, offset_s};
		m_windowHasPoint = true;
	}
	else if (offset_s < m_windowMin.offset_s)
	{
		m_windowMin = Point{sensor_s, offset_s};
	}

	const double stamp_s = sensor_s + predictOffset(sensor_s) + slewCorrection(sensor_s);
	return m_hostRefNs + static_cast<int64_t>(std::llround(stamp_s * 1e9));
}

void HostClockSync::start(uint32_t sample_time_fine, int64_t host_ns)
{
	m_started = true;
	m_prevTicks = sample_time_fine;
	m_ticks = 0;
	m_hostRefNs = host_ns;
}

void HostClockSync::closeWindow(double now_sensor_s)
{
	m_minima.push_back(m_windowMin);
	while (m_minima.size() > m_config.history)
		m_minima.pop_front();
	m_windowHasPoint = false;
	if (m_delayCount > 0)
		m_meanDelay = m_delaySum / static_cast<double>(m_delayCount);
	m_delaySum = 0.0;
	m_delayCount = 0;
	// Blend the new line in: the stamps continue from where the old line (and
	// any blend still running) had them, and reach the new line over slew_s.
	const double before = predictOffset(now_sensor_s) + slewCorrection(now_sensor_s);
	const bool had_fit = m_windows > 0;
	++m_windows;
	fit();
	if (had_fit)
	{
		m_slewStart = now_sensor_s;
		m_slewAmount = before - predictOffset(now_sensor_s);
	}
}

void HostClockSync::fit()
{
	// Until the minima span enough time the slope is mostly noise: estimate
	// the offset only (130 ppm over 10 s is 1.3 ms).
	const Point &first = m_minima.front();
	const Point &last = m_minima.back();
	double lowest = first.offset_s;
	for (const Point &p : m_minima)
		lowest = std::min(lowest, p.offset_s);
	m_a = lowest;
	m_b = 0.0;
	if (last.sensor_s - first.sensor_s < m_config.min_drift_span_s)
		return;

	// Lower convex hull of the minima (they are sorted by sensor time).
	std::vector<Point> hull;
	hull.reserve(m_minima.size());
	for (const Point &p : m_minima)
	{
		while (hull.size() >= 2)
		{
			const Point &o = hull[hull.size() - 2];
			const Point &a = hull.back();
			const double cross = (a.sensor_s - o.sensor_s) * (p.offset_s - o.offset_s) -
				(a.offset_s - o.offset_s) * (p.sensor_s - o.sensor_s);
			if (cross > 0.0)
				break;
			hull.pop_back();
		}
		hull.push_back(p);
	}

	// The line below all minima that minimizes the summed distance to them is
	// the hull edge spanning their mean sensor time. Only the fastest packets
	// define it, and it moves only when an even faster one arrives.
	double mean_s = 0.0;
	for (const Point &p : m_minima)
		mean_s += p.sensor_s;
	mean_s /= static_cast<double>(m_minima.size());
	std::size_t i = 0;
	while (i + 2 < hull.size() && hull[i + 1].sensor_s < mean_s)
		++i;
	const Point &p0 = hull[i];
	const Point &p1 = hull[std::min(i + 1, hull.size() - 1)];
	if (p1.sensor_s <= p0.sensor_s)
		return;
	const double b = (p1.offset_s - p0.offset_s) / (p1.sensor_s - p0.sensor_s);
	if (std::abs(b) > m_config.max_drift)
		return;  // not a crystal drift: keep the offset-only fit
	m_b = b;
	m_a = p0.offset_s - b * p0.sensor_s;
}

double HostClockSync::predictOffset(double sensor_s) const
{
	if (m_windows == 0)
	{
		// Nothing fitted yet: the fastest packet seen so far is the best guess.
		return m_windowHasPoint ? m_windowMin.offset_s : 0.0;
	}
	return m_a + m_b * sensor_s;
}

double HostClockSync::slewCorrection(double sensor_s) const
{
	if (m_slewAmount == 0.0)
		return 0.0;
	const double left = 1.0 - (sensor_s - m_slewStart) / m_config.slew_s;
	return left > 0.0 ? m_slewAmount * std::min(left, 1.0) : 0.0;
}

bool HostClockSync::locked() const
{
	return m_windows >= 2;
}

double HostClockSync::drift() const
{
	return m_b;
}

double HostClockSync::meanDelay() const
{
	return m_meanDelay;
}

std::size_t HostClockSync::resets() const
{
	return m_resets;
}

std::size_t HostClockSync::windows() const
{
	return m_windows;
}
