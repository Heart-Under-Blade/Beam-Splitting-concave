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
