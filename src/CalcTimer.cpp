#include "CalcTimer.h"

#include <limits.h>

#define MSEC_PER_SEC 1000
#define SEC_PER_MIN 60
#define MIN_PER_HOUR 60
#define HOUR_PER_DAY 24

using namespace std::chrono;

CalcTimer::CalcTimer()
{
	Reset();
}

void CalcTimer::Start()
{
	m_startPoint = system_clock::now();
}

void CalcTimer::Stop()
{
	m_nowPoint = system_clock::now();
}

long long CalcTimer::Duration()
{
	auto dur = m_nowPoint - m_startPoint;
	return duration_cast<milliseconds>(dur).count();
}

void CalcTimer::Left(const long long &ms)
{
	if (ms >= m_lastTimeLeft)
	{
		return; // time left must not increase
	}

	m_lastTimeLeft = ms;
	seconds = ms/MSEC_PER_SEC;

	auto SetUnitValue = [](int &unitVal, int &lessUnitVal, int lim)
	{
		if (lessUnitVal > lim)
		{
			unitVal = lessUnitVal/lim;
			lessUnitVal %= lim;
		}
		else
		{
			unitVal = 0;
		}
	};

	SetUnitValue(minutes, seconds, SEC_PER_MIN);
	SetUnitValue(hours, minutes, MIN_PER_HOUR);
	SetUnitValue(days, hours, HOUR_PER_DAY);
}

const time_t &CalcTimer::End(const long long &ms)
{
	if (ms <= m_lastTimeLeft) // end time must not increase
	{
		auto now = m_nowPoint;
		now += milliseconds(m_lastTimeLeft);
		m_lastEndTime = system_clock::to_time_t(now);
	}

	return m_lastEndTime;
}

void CalcTimer::Reset()
{
	days = 0;
	hours = 0;
	minutes = 0;
	seconds = 0;
	m_lastTimeLeft = LONG_MAX;
	m_lastEndTime = 0;
}

std::string CalcTimer::ToString()
{
	std::string strTime;

	if (days != 0) // REF: merge
	{
		strTime.append(std::to_string(days));
		strTime.append("d ");
	}

	if (hours != 0)
	{
		strTime.append(std::to_string(hours));
		strTime.append("h ");
	}

	if (minutes != 0)
	{
		strTime.append(std::to_string(minutes));
		strTime.append("m ");
	}

	if (seconds != 0)
	{
		strTime.append(std::to_string(seconds));
		strTime.append("s ");
	}

	return strTime;
}
