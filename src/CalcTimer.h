#pragma once

#include <chrono>
#include <string>

class CalcTimer
{
public:
	CalcTimer();

	void Start();
	void Stop();
	long long Duration();
	void Left(const long long &ms);
	const std::time_t &End(const long long &ms);
	void Reset();
	std::string ToString();

public:
	int days;
	int hours;
	int minutes;
	int seconds;

private:
	std::chrono::time_point<std::chrono::system_clock> m_startPoint;
	std::chrono::time_point<std::chrono::system_clock> m_nowPoint;
	long long m_lastTimeLeft;
	std::time_t m_lastEndTime;
};
