#pragma once

#include "Scattering.h"

class ScatteringConvex : public Scattering
{
public:
	ScatteringConvex(Particle *particle, Light *incidentLight,
					 bool isOpticalPath, int nActs);

	void ScatterLight(std::vector<Beam> &outBeams) override;
	void ScatterLight(const std::vector<std::vector<int>> &,
					  std::vector<Beam> &) override; ///> for predefined trajectories

protected:
	void TraceInternalBeams(std::vector<Beam> &outBeams);
	void SplitSecondaryBeams(Beam &incidentBeam, Facet *facet,
							 std::vector<Beam> &outBeams);
};
