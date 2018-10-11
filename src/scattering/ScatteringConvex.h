#pragma once

#include "Scattering.h"

class ScatteringConvex : public Scattering
{
public:
	ScatteringConvex(Particle *particle, Light *incidentLight,
					 bool isOpticalPath, int nActs);
protected:
	void SplitBeams(std::vector<Beam> &scatteredBeams) override;
	void SplitLightToBeams(std::vector<Beam> &scatteredBeams) override;

	void SplitBeamByFacets(Beam &beam, std::vector<Beam> &scatteredBeams);
};
