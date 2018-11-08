#pragma once

#include "RegularIncidence.h"

class NormalIncidence : public RegularIncidence
{
public:
	void ComputeDirections(Beam &beam, SplittedBeams<Beam> &beams,
						   Splitting &splitter) const;
	void ComputeJonesMatrices(Beam &beam, SplittedBeams<Beam> &beams,
							  Splitting &splitter) const;
};
