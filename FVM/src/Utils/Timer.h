#pragma once

#include <chrono>


class Timer
{
public:

    void start();
    void stop();

    // get elapsed time in seconds
    double getElapsedTime() const;

private:

    using TimerType     = std::chrono::steady_clock;
    using DurationType  = TimerType::duration;
    using TimePointType = std::chrono::time_point<TimerType>;

    DurationType  m_elapsed = DurationType::zero();
    TimePointType m_start_time;
};
