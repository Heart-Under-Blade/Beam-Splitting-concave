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

std::string CalcTimer::Start()
{
	m_startPoint = system_clock::now();
	time_t t = system_clock::to_time_t(m_startPoint);
	m_startTime = ctime(&t);
	return m_startTime;
}

std::string CalcTimer::Stop()
{
	m_nowPoint = system_clock::now();
	time_t t = system_clock::to_time_t(m_nowPoint);
	m_stopTime = ctime(&t);
	return m_stopTime;
}

std::string CalcTimer::Elapsed()
{
	m_nowPoint = system_clock::now();
	auto dur = m_nowPoint - m_startPoint;
	long long sec = duration_cast<seconds>(dur).count();

	m_seconds = sec;

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

	SetUnitValue(m_minutes, m_seconds, SEC_PER_MIN);
	SetUnitValue(m_hours, m_minutes, MIN_PER_HOUR);
	SetUnitValue(m_days, m_hours, HOUR_PER_DAY);

	return ToString();
}

std::string CalcTimer::GetStartTime() const
{
	return m_startTime;
}

std::string CalcTimer::GetStopTime() const
{
	return m_stopTime;
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
	m_seconds = ms/MSEC_PER_SEC;

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

	SetUnitValue(m_minutes, m_seconds, SEC_PER_MIN);
	SetUnitValue(m_hours, m_minutes, MIN_PER_HOUR);
	SetUnitValue(m_days, m_hours, HOUR_PER_DAY);
}

void CalcTimer::Reset()
{
	m_days = 0;
	m_hours = 0;
	m_minutes = 0;
	m_seconds = 0;
	m_lastTimeLeft = LONG_MAX;
	m_lastEndTime = 0;
}

std::string CalcTimer::ToString()
{
	std::string strTime;

	if (m_days != 0) // REF: merge
	{
		strTime.append(std::to_string(m_days));
		strTime.append("d ");
	}

	if (m_hours != 0)
	{
		strTime.append(std::to_string(m_hours));
		strTime.append("h ");
	}

	if (m_minutes != 0)
	{
		strTime.append(std::to_string(m_minutes));
		strTime.append("m ");
	}

	if (m_seconds != 0)
	{
		strTime.append(std::to_string(m_seconds));
		strTime.append("s ");
	}

	return strTime;
}
