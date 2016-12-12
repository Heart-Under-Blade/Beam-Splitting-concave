#pragma once

#include <fstream>
static std::ofstream logfile("log.txt", std::ios::out);
static int assertNum = 0;

#define LOG_ASSERT(expr) \
if (!(expr)) \
{ \
	logfile << "ASSERT: \"" << #expr << "\", file: " << __FILE__ << ", line: " << __LINE__ << std::endl; \
	logfile.flush(); \
	++assertNum; \
}
