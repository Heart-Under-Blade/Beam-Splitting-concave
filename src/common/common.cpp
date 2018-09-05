#include "global.h"
#include "macro.h"
#include <string>

#include <iostream>
using namespace std;

#ifdef _WIN32
#include <windows.h>
#endif

string CreateUniqueFileName(const string &filename)
{
	string name = filename + ".dat";

	for (int i = 1; ifstream(name) != NULL; ++i)
	{
		name = filename + '(' + to_string(i) + ')' + ".dat";
	}

	return name;
}

string CreateFolder(string &name)
{
#ifdef _WIN32
	char curDir[MAX_PATH] = ""; // current directory
	char newDir[MAX_PATH]; // created directory

	// get the current directory, and store it
	if (!GetCurrentDirectoryA(MAX_PATH, curDir))
	{
		cerr << "Error getting current directory: #" << GetLastError();
	}

	strcat(curDir, "\\");
	strcpy(newDir, curDir);
	strcat(newDir, name.c_str());

	char dirN[MAX_PATH];
	strcpy(dirN, newDir);
	string num = "";

	for (int i = 1; !CreateDirectoryA(dirN, NULL); ++i)
	{
		num = "(" + to_string(i) + ")";
		string dirName = newDir + num;
		strcpy(dirN, dirName.c_str());
	}

	name += num;
#else
	cr_dir = name;
#endif
	return curDir;
}

string CreateDir(const string &name)
{
	string dirName;
#ifdef _WIN32
	char cdir[MAX_PATH];
	cdir[0] = '\0';
	strcat(cdir, name.c_str());

	char dirN[MAX_PATH];
	strcpy(dirN, cdir);

	for (int i = 1; !CreateDirectoryA(cdir, NULL); ++i)
	{
		string dirName = dirN;
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
