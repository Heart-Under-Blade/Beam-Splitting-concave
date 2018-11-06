#pragma once

#include "Beam.h"
#include "PathedBeam.h"
#include "Splitting.h"

class Incidence
{
public:
	virtual void ComputeDirections(Beam &/*beam*/, SplittedBeams<Beam> &/*beams*/,
								   Splitting &/*splitter*/) const {};
	virtual void ComputeJonesMatrices(Beam &/*beam*/, SplittedBeams<Beam> &/*beams*/,
									  Splitting &/*splitter*/) const {};

	virtual void ComputeOpticalPaths(const PathedBeam &beam, SplittedBeams<PathedBeam> &beams,
									 Splitting &splitter) const;
};
