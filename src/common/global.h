#pragma once

#define M_PI	3.14159265358979323846
#define M_2PI	6.283185307179586476925286766559

// debug
#define _TRACK_ALLOW // tracks of beams writes to beams
//#define _CALC_AREA_CONTIBUTION_ONLY
#define _CHECK_ENERGY_BALANCE // energy conversation of input and output beams outputs in file

#include <string>

void EraseConsoleLine(int lenght);

std::string CreateFolder(std::string &name); // returned path to folder and created folder name
std::string CreateDir(const std::string &name);
std::string CreateUniqueFileName(const std::string &filename);
