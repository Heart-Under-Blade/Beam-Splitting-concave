#include "global.h"
#include "macro.h"
#include <string>

#include <iostream>
using namespace std;

#ifdef _WIN32
#include <windows.h>
#endif

std::string CreateDir(const std::string &name)
{
	std::string dirName;

#ifdef _WIN32
	char dir[260] = "";
	size_t bufferSize = MAX_PATH;
	char cdir[MAX_PATH]; // store the current directory

	// get the current directory, and store it
	if (!GetCurrentDirectoryA(bufferSize, cdir))
	{
		cerr << "Error getting current directory: #" << GetLastError();
	}

	strcat(cdir, "\\");
	strcat(cdir, name.c_str());

	char dirN[MAX_PATH];
	strcpy(dirN, cdir);

	for (int i = 1; !CreateDirectoryA(cdir, NULL); ++i)
	{
		std::string dirName = dirN;
		dirName += "(" + to_string(i) + ")";
		strcpy(cdir, dirName.c_str());
	}

	strcat(cdir, "\\");
	strcat(dir, cdir);
	dirName = dir;
#else
	dirName = name;
#endif
	return dirName;
}

std::string CreateDir2(const std::string &name)
{
	std::string dirName;
#ifdef _WIN32
	char cdir[MAX_PATH];
	cdir[0] = '\0';
	strcat(cdir, name.c_str());

	char dirN[MAX_PATH];
	strcpy(dirN, cdir);

	for (int i = 1; !CreateDirectoryA(cdir, NULL); ++i)
	{
		std::string dirName = dirN;
		dirName += "(" + to_string(i) + ")";
		strcpy(cdir, dirName.c_str());
	}

	strcat(cdir, "\\");
	dirName = cdir;
#else
	dirName = name;
#endif
	return dirName;
}

void EraseConsoleLine(int lenght)
{
	cout << '\r';

	for (int i = 0; i < lenght; ++i)
	{
		cout << ' ';
	}

	cout << '\r';
}

double DegToRad(double deg)
{
	return (deg*M_PI)/180;
}

double RadToDeg(double rad)
{
	return (rad*180)/M_PI;
}
