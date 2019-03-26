#pragma once

#include "ArgPP.h"
#include "geometry_lib.h"

struct AngleRange
{
	Orientation min;
	Orientation max;
	int numberZenith;
	int numberAzimuth;
	double normZenith;
	double stepZenith;
	double normAzimuth;
	double stepAzimuth;

    AngleRange() {}

	void Setup()
	{
		normZenith = max.zenith - min.zenith;
		stepZenith = normZenith/numberZenith;

		normAzimuth = max.azimuth - min.azimuth;
		stepAzimuth = normAzimuth/numberAzimuth;
    }
};

class ArgumentParser
{
public:
    ArgumentParser(ArgPP *parser);

    void SetArgRules();
    AngleRange GetRange(const Orientation &symmetry);

private:
    ArgPP *m_parser;

    void GetRangeParams(const Orientation &symmetry, const std::string &key,
                        double &min, double &max);
};
