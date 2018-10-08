#pragma once

#include "Scattering.h"

class ScatteringConvex : public Scattering
{
public:
	ScatteringConvex(Particle *particle, Light *incidentLight,
					 bool isOpticalPath, int nActs);

	void ScatterLight(std::vector<Beam> &scatteredBeams) override;
	void ScatterLight(const std::vector<std::vector<int>> &,
					  std::vector<Beam> &) override; ///> for predefined trajectories

protected:
	void ScatterBeams(std::vector<Beam> &outBeams);
	void SplitBeamsByFacet(Beam &beam, Facet *facet,
						   std::vector<Beam> &outBeams);
	void SplitLightToBeams(std::vector<Beam> &scatteredBeams);
};
