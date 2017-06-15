#pragma once

#include <fstream>
static std::ofstream logfile("log.txt", std::ios::out);

#define LOG_ASSERT(expr) \
if (!(expr)) \
{ \
	logfile << "ASSERT: \"" << #expr << "\", file: " << __FILE__ << ", line: " << __LINE__; \
	throw false; \
}

void OutputState(int i, int j);
