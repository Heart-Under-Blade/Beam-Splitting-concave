#include "global.h"
#include "macro.h"
#include <string>
#include <iostream>

#ifdef _WIN32
#include <windows.h>
#endif

void OutputState(int i, int j)
{
	logfile << "i: " << i << "; j: " << j << std::endl;
	logfile.flush();
}

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

void CreateDir(const char *name, char *dir)
{
#ifdef _WIN32
	size_t bufferSize = MAX_PATH;
	char cdir[MAX_PATH]; // store the current directory

	// get the current directory, and store it
	if (!GetCurrentDirectoryA(bufferSize, cdir))
	{
		std::cerr << "Error getting current directory: #" << GetLastError();
	}

	strcat(cdir, "\\");
	strcat(cdir, name);

	char dirN[MAX_PATH];
	strcpy(dirN, cdir);

	for (int i = 1; !CreateDirectoryA(cdir, NULL); ++i)
	{
		std::string dirName = dirN;
		dirName += "(" + std::to_string(i) + ")";
		strcpy(cdir, dirName.c_str());
	}

//	if (!CreateDirectoryA(cdir, NULL))
//	{
//		std::cerr << "Error creating directory: #" << GetLastError();
//	}

	strcat(cdir, "\\");
	strcat(dir, cdir);
#else
	strcpy(dir, name);
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
	return (rad*180)/M_PI;
}
