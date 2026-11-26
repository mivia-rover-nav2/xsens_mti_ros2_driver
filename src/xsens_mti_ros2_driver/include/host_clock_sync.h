#ifndef HOST_CLOCK_SYNC_H
#define HOST_CLOCK_SYNC_H

#include <cstddef>
#include <cstdint>
#include <deque>

/*
 * Maps the MTi's SampleTimeFine counter onto the host clock.
 *
 * Each packet gives two times for the same sample:
 * - SampleTimeFine: the sampling instant on the MTi's own clock (10 kHz ticks).
 *   Its spacing is exact, but the clock has an unknown offset from the host
 *   and runs at a slightly different rate (drift, ~100 ppm on an MTi-8);
 * - the host clock when the packet reached the driver: the right clock, but
 *   late by the transport delay, which varies by tens of milliseconds because
 *   the data arrives over USB in bursts.
 *
 * The delay is never negative, and the packets that arrived fastest are the
 * closest to the true sampling instant. For each window of sensor time the
 * class keeps the packet with the smallest (host - sensor) difference, and fits
 * the line that runs below the last windows' minima and touches the lowest ones
 * (an edge of their lower convex hull):
 *
 *     host = sensor + a + b * sensor
 *
 * where a absorbs the offset and b the drift (estimated only once the minima
 * span 10 s; before that, the offset alone). Stamps then keep the sensor's
 * exact spacing while living on the host clock; when the line changes, the
 * difference is blended in over 2 s instead of stepping the stamps. A step of either clock (NTP
 * correcting the host, the MTi restarting) is detected and the fit restarts.
 *
 * Pure C++ (no ROS) so it can be unit tested with simulated clocks.
 */
class HostClockSync
{
public:
	struct Config
	{
		double window_s = 1.0;           // length of one lower-envelope window, sensor time
		std::size_t history = 60;        // number of windows the line is fitted on
		double step_threshold_s = 0.1;   // a residual beyond this means a clock stepped...
		double step_confirm_s = 0.25;    // ...once it lasts this long above the line
		double min_drift_span_s = 10.0;  // fit the drift only once the minima span this long
		double max_drift = 1e-3;         // larger slopes are not a crystal drift
		double slew_s = 2.0;             // a refit is blended into the stamps over this long
	};

	HostClockSync();
	explicit HostClockSync(const Config &config);

	// sample_time_fine: raw SampleTimeFine (10 kHz ticks, wraps at 2^32).
	// host_ns: host clock when the packet arrived, nanoseconds.
	// Returns the sampling instant on the host clock, nanoseconds.
	int64_t update(uint32_t sample_time_fine, int64_t host_ns);

	void reset();

	// True once the line is fitted on at least two windows.
	bool locked() const;
	// Drift of the host clock relative to the sensor clock (1e-6 = 1 ppm).
	double drift() const;
	// Mean transport delay above the fastest packets in the last closed window,
	// seconds (how bursty the arrival is).
	double meanDelay() const;
	// Number of restarts after a clock step.
	std::size_t resets() const;
	// Number of windows closed since the last restart.
	std::size_t windows() const;

private:
	struct Point
	{
		double sensor_s;  // sensor time of the fastest packet in the window
		double offset_s;  // host - sensor for that packet
	};

	void clear();
	void start(uint32_t sample_time_fine, int64_t host_ns);
	void closeWindow(double now_sensor_s);
	void fit();
	double predictOffset(double sensor_s) const;
	double slewCorrection(double sensor_s) const;

	Config m_config;
	bool m_started;
	uint32_t m_prevTicks;
	int64_t m_ticks;        // unwrapped ticks since the first sample
	int64_t m_hostRefNs;    // host time of the first sample

	long m_window;          // index of the window being filled
	bool m_windowHasPoint;
	Point m_windowMin;
	bool m_aboveLine;       // in a run of samples far above the fitted line
	double m_aboveSince;    // sensor time the run started
	double m_delaySum;      // residuals above the line in the current window
	std::size_t m_delayCount;
	double m_meanDelay;
	std::deque<Point> m_minima;

	double m_a;
	double m_b;
	double m_slewStart;     // sensor time of the last refit
	double m_slewAmount;    // stamp correction at m_slewStart, fading to 0 over slew_s
	std::size_t m_resets;
	std::size_t m_windows;
};

#endif // HOST_CLOCK_SYNC_H
