#pragma once

#include "HandlerGO.h"

class HandlerTracksGO : public HandlerGO
{
public:
	HandlerTracksGO(Particle *particle, Light *incidentLight, float wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
	void WriteMatricesToFile(std::string &destName) override;
};
