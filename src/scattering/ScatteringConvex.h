#pragma once

#include "Scattering.h"

class ScatteringConvex : public Scattering
{
public:
	ScatteringConvex(Particle *particle, Light *incidentLight,
					 bool isOpticalPath, int nActs);

	void ScatterLight(double beta, double gamma, std::vector<Beam> &outBeams) override;
	void ScatterLight(double, double, const std::vector<std::vector<int>> &,
					  std::vector<Beam> &) override; ///> for predefined trajectories

protected:
	void TraceInternalBeams(std::vector<Beam> &outBeams);

	bool SplitSecondaryBeams(Beam &incidentBeam, int facetID,
							 Beam &inBeam, std::vector<Beam> &outBeams);
};
