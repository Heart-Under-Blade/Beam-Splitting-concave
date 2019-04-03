#pragma once

#include <string>
#include <map>
#include <vector>

class ScatteringFiles
{
public:
	ScatteringFiles(const std::string &dir, const std::string &folder,
					const std::string &tableHead = "");

	void CreateMainFile(const std::string &subdir, const std::string &name);
	void CreateGroupFile(const std::string &subdir, const std::string &name);

	std::ofstream *GetMainFile(const std::string &name);
	std::ofstream *GetGroupFile(int i);

	~ScatteringFiles();

private:
	std::map<std::string, std::ofstream*> m_mainFiles;
	std::vector<std::ofstream*> m_groupFiles;

	std::string m_dir;
	std::string m_folder;
	std::string m_fullDir;
	std::string m_tableHead;

	std::ofstream *CreateFile(const std::string &name,
							  const std::string &subfolder = "");
};
