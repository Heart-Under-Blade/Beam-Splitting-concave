#pragma once

#include "Tracer.h"

class TracerGO : public Tracer
{
public:
	TracerGO(Particle *particle, Scattering *scattering,
			 const std::string &resultFileName);

	void TraceRandom(const OrientationRange &range) override;

protected:
	double CalcNorm(long long orNum);
	void OutputSummary(int orNumber, double D_tot, double NRM, CalcTimer &timer);
};
