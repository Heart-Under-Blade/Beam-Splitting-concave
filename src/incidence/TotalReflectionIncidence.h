#pragma once

#include "Incidence.h"

class TotalReflectionIncidence : public Incidence
{
public:
	TotalReflectionIncidence();
	void ComputeDirections(Beam &beam, BeamPair<Beam> &beams) const override;
	void ComputeJonesMatrices(Beam &parentBeam, BeamPair<Beam> &beams) const override;

	void ComputeOpticalPaths(const PathedBeam &parentBeam,
							 BeamPair<PathedBeam> &beams) const override;
};
