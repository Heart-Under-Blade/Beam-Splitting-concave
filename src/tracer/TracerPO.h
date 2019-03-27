#pragma once

#include "Tracer.h"

class TracerPO : public LightTracer
{
public:
	TracerPO(Particle *particle, Scattering *scattering,
			 const std::string &resultFileName);

	void TraceRandom(const AngleRange &zenithRange,
					 const AngleRange &azimuthRange);
};
