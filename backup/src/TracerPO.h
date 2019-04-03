#pragma once

#include "Tracer.h"

class TracerPO : public Tracer
{
public:
	TracerPO(Particle *particle, int nActs, const std::string &resultFileName);
	void TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange);
	void TraceFixed(const double &beta, const double &gamma);
};
