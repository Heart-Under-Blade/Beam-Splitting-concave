#pragma once

#include "Incidence.h"

class CompleteReflectionIncidence : public Incidence
{
public:
	void ComputeDirections(Beam &beam, SplittedBeams<Beam> &beams,
						   Splitting &splitter) const override;
	void ComputeJonesMatrices(Beam &parentBeam, SplittedBeams<Beam> &beams,
							  Splitting &splitter) const override;

	void ComputeOpticalPaths(const PathedBeam &incidentBeam,
							 SplittedBeams<PathedBeam> &beams,
							 Splitting &splitter) const override;
};
