#include "global.h"
#include <string>
#include <iostream>

#ifdef _WIN32
#include <windows.h>
#endif

void Dellines(int count)
{
#ifdef _WIN32
	std::string mask;

	for (int i = 0; i < count*80; ++i)
	{
		mask.append(" ");
	}

	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO ci;
	GetConsoleScreenBufferInfo(hConsole, &ci);

	ci.dwCursorPosition.X = 0;
	ci.dwCursorPosition.Y -= count-1;
	SetConsoleCursorPosition(hConsole, ci.dwCursorPosition);
	std::cout<< mask;
	SetConsoleCursorPosition(hConsole, ci.dwCursorPosition);
#else
	++count;
	--count;
#endif
}

void EraseConsoleLine(int lenght)
{
	std::cout << '\r';

	for (int i = 0; i < lenght; ++i)
	{
		std::cout << ' ';
	}

	std::cout << '\r';
}

double DegToRad(double deg)
{
	return (deg*M_PI)/180;
}

double RadToDeg(double rad)
{
	rad = (rad*180);
	rad = rad/M_PI;
	return rad;
}
