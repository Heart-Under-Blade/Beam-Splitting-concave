#pragma once

#include <chrono>
#include <string>

class CalcTimer
{
public:
	CalcTimer();

	time_t Start();
	time_t Stop();
	long long Duration();
	void Left(const long long &ms);
	const std::time_t &End(const long long &ms);
	time_t Begin() const;
	void Reset();
	std::string ToString();
	std::string Elapsed();
	long long SecondsElapsed();

private:
	std::chrono::time_point<std::chrono::system_clock> m_startPoint;
	std::chrono::time_point<std::chrono::system_clock> m_nowPoint;
	long long m_lastTimeLeft;
	std::time_t m_lastEndTime;

	int m_days;
	int m_hours;
	int m_minutes;
	int m_seconds;
};
