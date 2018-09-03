#include "ScatteringFiles.h"

#include <fstream>
#include <iomanip>
#include <assert.h>
#include <limits>

#define RES_EXT ".dat"

using namespace std;

ScatteringFiles::ScatteringFiles(const string &dir, const string &folder,
								 const string &tableHead)
	: m_dir(dir), m_folder(folder), m_tableHead(tableHead)
{
	m_fullDir = m_dir + m_folder + '\\';
}

void ScatteringFiles::CreateMainFile(const string &subdir, const string &name)
{
	ofstream *file = CreateFile(name, subdir);
	pair<string, ofstream*> fileItem(name, file);
	m_mainFiles.insert(fileItem);
#ifdef _DEBUG // DEB
//	auto it = m_mainFiles.find(name);
//	bool d = it != m_mainFiles.end();
//	ofstream *dd = m_mainFiles[name];
//	*(it->second) << (2232);
//	int ff = 0;
#endif
}

void ScatteringFiles::CreateGroupFile(const string &subdir, const string &name)
{
	ofstream *file = CreateFile(name, subdir);
	m_groupFiles.push_back(file);
}

ofstream *ScatteringFiles::GetMainFile(const string &name)
{
	return m_mainFiles[name];
}

ofstream *ScatteringFiles::GetGroupFile(int i)
{
	return m_groupFiles[i];
}

ofstream *ScatteringFiles::CreateFile(const string &name, const string &subfolder)
{
	string filename = m_fullDir + subfolder + '\\' + name + "__" + m_folder + RES_EXT;
	ofstream *file = new ofstream(filename, ios::out);
	assert(file->is_open());
	(*file) << setprecision(numeric_limits<long double>::digits10 + 1);
	(*file) << m_tableHead;
	return file;
}

ScatteringFiles::~ScatteringFiles()
{
	for (const auto &p : m_mainFiles)
	{
		p.second->close();
		delete p.second;
	}

	for (const auto &p : m_groupFiles)
	{
		p->close();
		delete p;
	}
}
