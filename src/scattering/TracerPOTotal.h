#pragma once

#include "TracerPO.h"

class TracerPOTotal : public TracerPO
{
public:
	TracerPOTotal(Particle *particle, Scattering *scattering,
				  const std::string &resultFileName);
	void TraceRandom(const OrientationRange &range) override;
};
