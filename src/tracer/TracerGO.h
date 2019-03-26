#pragma once

#include "Tracer.h"

class TracerGO : public LightTracer
{
public:
	TracerGO(Particle *particle, Scattering *scattering,
			 const std::string &resultFileName);

	void TraceRandom(const AngleRange &zenithRange,
					 const AngleRange &azimuthRange);

protected:
	double CalcNorm(long long orNum);
	void OutputSummary(int orNumber, double D_tot, double NRM, CalcTimer &timer);
};
