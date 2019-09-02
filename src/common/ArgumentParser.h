#pragma once

#include "ArgPP.h"
#include "geometry_lib.h"

class ArgumentParser
{
public:
    ArgumentParser(ArgPP *parser);

    void SetArgRules();
	OrientationRange GetRange(const Orientation &symmetry);

private:
    ArgPP *m_parser;

    void GetRangeParams(const Orientation &symmetry, const std::string &key,
                        double &min, double &max);
};
