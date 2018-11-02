#pragma once

#include "RegularIncidence.h"

class Beam;
class Splitting;

class CompleteReflectionIncidence : public RegularIncidence
{
public:
	static void ComputeDirections(Beam &beam, Splitting &splitter);
	static void ComputeJonesMatrices(Beam &parentBeam, SplittedBeams<Beam> &beams,
									  Splitting &splitter);

	static void ComputeOpticalPaths(const PathedBeam &incidentBeam,
									SplittedBeams<PathedBeam> &beams,
									Splitting &splitter);
};
