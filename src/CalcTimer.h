#pragma once

#include <chrono>
#include <string>

class CalcTimer
{
public:
	CalcTimer();

	std::string Start();
	std::string Stop();
	long long Duration();
	void Left(const long long &ms);
	void Reset();
	std::string ToString();
	std::string Elapsed();

	std::string GetStartTime() const;
	std::string GetStopTime() const;

private:
	std::chrono::time_point<std::chrono::system_clock> m_startPoint;
	std::chrono::time_point<std::chrono::system_clock> m_nowPoint;

	std::string m_startTime;
	std::string m_stopTime;

	long long m_lastTimeLeft;
	std::time_t m_lastEndTime;

	int m_days;
	int m_hours;
	int m_minutes;
	int m_seconds;
};
