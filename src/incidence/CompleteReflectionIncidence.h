#pragma once

#include "Incidence.h"

class CompleteReflectionIncidence : public Incidence
{
public:
	CompleteReflectionIncidence();
	void ComputeDirections(Beam &beam, BeamPair<Beam> &beams) const override;
	void ComputeJonesMatrices(Beam &parentBeam, BeamPair<Beam> &beams) const override;

	void ComputeOpticalPaths(const PathedBeam &incidentBeam,
							 BeamPair<PathedBeam> &beams) const override;
};
