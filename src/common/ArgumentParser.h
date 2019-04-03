#pragma once

#include "ArgPP.h"
#include "geometry_lib.h"

struct AngleRange
{
	double min;
	double max;
	int number;
	double norm;
	double step;

    AngleRange() {}

	void Setup()
	{
		norm = max - min;
		step = norm/number;
    }
};

class ArgumentParser
{
public:
    ArgumentParser(ArgPP *parser);

    void SetArgRules();
	AngleRange GetRange(const std::string &key,
						const Orientation &symmetry);
private:
    ArgPP *m_parser;
    void GetRangeParams(const Orientation &symmetry, const std::string &key,
                        double &min, double &max);
};
