// Simulated MTi stream: 100 Hz samples, the sensor clock 130 ppm slower than
// the host and started at an arbitrary counter value, packets reaching the
// host in USB bursts every 40 ms. HostClockSync must put the stamps on the
// host clock with a constant delay (the arrival time jitters by 40 ms), keep
// the 10 ms spacing, and restart after clock steps.
//
// The constant part cannot be removed from arrival times alone: it is the
// delay of the fastest packet, here up to one sample period because a burst
// boundary falls at most 10 ms after a sample.
#include "host_clock_sync.h"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

static void check(bool ok, const std::string &message)
{
	if (!ok)
		throw std::runtime_error(message);
}

namespace
{
constexpr double kRate = 100.0;
constexpr double kDrift = 130e-6;          // sensor clock slower than the host
constexpr int64_t kHostEpochNs = 1790000000LL * 1000000000LL;
constexpr double kBurst = 0.040;           // USB hands data over every 40 ms
constexpr double kMinDelay = 0.001;

struct Stream
{
	std::mt19937 rng{42};
	std::uniform_real_distribution<double> jitter{0.0, 0.0005};
	uint32_t ticks0 = 0xFFFFFFFFu - 200000u;  // the counter wraps after ~20 s
	double host_step = 0.0;                    // host clock step applied so far, s
	double tick_restart_at = -1.0;             // true time the MTi restarted

	uint32_t ticks(double t) const
	{
		if (tick_restart_at >= 0.0 && t >= tick_restart_at)
			return static_cast<uint32_t>(std::floor((t - tick_restart_at) * 1e4 * (1.0 - kDrift)));
		return ticks0 + static_cast<uint32_t>(std::floor(t * 1e4 * (1.0 - kDrift)));
	}
	// Host time at which the packet sampled at true time t reaches the driver.
	double arrival(double t)
	{
		return std::ceil((t + kMinDelay) / kBurst) * kBurst + jitter(rng);
	}
};

int64_t ns(double s)
{
	return kHostEpochNs + static_cast<int64_t>(std::llround(s * 1e9));
}
}

int main()
{
	HostClockSync sync;
	Stream stream;

	double min_err = 1e9, max_err = -1e9, max_naive_err = 0.0, max_spacing_err = 0.0, max_detection = 0.0;
	double prev_stamp = 0.0;
	bool have_prev = false;
	const double settle = 15.0;  // warm-up (minima spanning 10 s, offset only) plus the 2 s blend
	double last_event_t = 0.0;
	bool waiting_detection = false;
	std::size_t seen_resets = sync.resets();

	const int n = static_cast<int>(140 * kRate);
	for (int k = 0; k < n; ++k)
	{
		const double t = k / kRate;
		bool event = true;
		if (k == static_cast<int>(60 * kRate))
			stream.host_step = 2.0;    // NTP steps the host clock forward
		else if (k == static_cast<int>(80 * kRate))
			stream.host_step = 2.2;    // a small forward step
		else if (k == static_cast<int>(100 * kRate))
			stream.host_step = 0.5;    // and back
		else if (k == static_cast<int>(120 * kRate))
			stream.tick_restart_at = t;  // the MTi restarts its counter
		else
			event = false;
		if (event)
		{
			last_event_t = t;
			waiting_detection = true;
			have_prev = false;
		}

		const double arrival = stream.arrival(t) + stream.host_step;
		const int64_t stamp_ns = sync.update(stream.ticks(t), ns(arrival));
		const double stamp = static_cast<double>(stamp_ns - kHostEpochNs) * 1e-9;
		const double truth = t + stream.host_step;

		if (sync.resets() != seen_resets)
		{
			seen_resets = sync.resets();
			check(waiting_detection, "restart without a clock step at t=" + std::to_string(t));
			max_detection = std::max(max_detection, t - last_event_t);
			waiting_detection = false;
		}
		if (t - last_event_t < settle)
			continue;
		check(!waiting_detection, "clock step not detected");

		const double err = stamp - truth;
		min_err = std::min(min_err, err);
		max_err = std::max(max_err, err);
		max_naive_err = std::max(max_naive_err, std::abs(arrival - truth));
		if (have_prev)
			max_spacing_err = std::max(max_spacing_err, std::abs((stamp - prev_stamp) - 1.0 / kRate));
		prev_stamp = stamp;
		have_prev = true;

		if (t > 20.0 && t < 60.0)
		{
			check(sync.locked(), "locked after the warm-up");
			check(std::abs(sync.drift() - kDrift) < 10e-6,
				"drift estimate " + std::to_string(sync.drift() * 1e6) + " ppm");
			check(sync.meanDelay() > 0.010 && sync.meanDelay() < 0.030,
				"mean delay " + std::to_string(sync.meanDelay()));
		}
	}

	std::cout << "stamp delay " << min_err * 1e3 << " to " << max_err * 1e3 << " ms (arrival time: up to "
		<< max_naive_err * 1e3 << " ms), max spacing error " << max_spacing_err * 1e6
		<< " us, resets " << sync.resets() << ", slowest step detection "
		<< max_detection * 1e3 << " ms" << std::endl;

	check(max_naive_err > 0.030, "the simulated arrival is bursty");
	check(max_err - min_err < 0.001, "the stamp delay is constant within 1 ms");
	check(min_err >= -0.001 && max_err < 1.0 / kRate + 0.002, "the constant delay is at most one sample period");
	check(max_spacing_err < 0.0002, "10 ms spacing kept within 0.2 ms");
	check(sync.resets() == 4, "one restart per step: " + std::to_string(sync.resets()));
	check(max_detection < 0.3, "steps detected within 0.3 s");
	std::cout << "host clock sync: all checks passed" << std::endl;
	return 0;
}
