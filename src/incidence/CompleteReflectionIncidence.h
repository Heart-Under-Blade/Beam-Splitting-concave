#pragma once

#include "Incidence.h"

class CompleteReflectionIncidence : public Incidence
{
public:
	CompleteReflectionIncidence(const complex &ri);
	void ComputeDirections(Beam &beam, SplittedBeams<Beam> &beams) const override;
	void ComputeJonesMatrices(Beam &parentBeam, SplittedBeams<Beam> &beams) const override;

	void ComputeOpticalPaths(const PathedBeam &incidentBeam,
							 SplittedBeams<PathedBeam> &beams) const override;

	double cosA2;
};
