#include "Timer.h"


void Timer::start()
{
    m_start_time = TimerType::now();
}


void Timer::stop()
{
    TimePointType end = TimerType::now();
    m_elapsed += end - m_start_time;
}


double Timer::getElapsedTime() const
{
    using namespace std::chrono;
    using TickType = microseconds;
    
    constexpr int64_t ticksPerSec = duration_cast<TickType>(seconds(1)).count();
    int64_t ticksElapsed = duration_cast<TickType>(m_elapsed).count();

    return  (double)ticksElapsed / ticksPerSec;
}
