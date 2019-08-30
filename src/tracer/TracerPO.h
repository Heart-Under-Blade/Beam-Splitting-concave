#pragma once

#include "Tracer.h"

class TracerPO : public Tracer
{
public:
	TracerPO(Particle *particle, Scattering *scattering,
			 const std::string &resultFileName);

	void TraceRandom(const OrientationRange &range) override;
};
