#pragma once

#include "HandlerGO.h"

class HandlerTotalGO : public HandlerGO
{
public:
	HandlerTotalGO(Particle *particle, Light *incidentLight, float wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
};
