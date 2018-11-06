#pragma once

#include "Incidence.h"

class RegularIncidence : public Incidence
{
public:
	void ComputeDirections(Beam &beam, SplittedBeams<Beam> &beams,
						   Splitting &splitter) const;
	void ComputeJonesMatrices(Beam &beam, SplittedBeams<Beam> &beams,
							  Splitting &splitter) const;
};
