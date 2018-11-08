#pragma once

#include "Incidence.h"
#include "PathedBeam.h"

class CompleteReflectionIncidence : public Incidence
{
public:
	void ComputeDirections(Beam &beam, SplittedBeams<Beam> &beams,
						   Splitting &splitter) const;

	void ComputeJonesMatrices(Beam &parentBeam, SplittedBeams<Beam> &beams,
							  Splitting &splitter) const;

	void ComputeOpticalPaths(const PathedBeam &incidentBeam,
							 SplittedBeams<PathedBeam> &beams,
							 Splitting &splitter) const;
};
