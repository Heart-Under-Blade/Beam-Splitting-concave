#include "global.h"
#include <windows.h>
#include <string>
#include <iostream>

void Dellines(int count)
{
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
}
