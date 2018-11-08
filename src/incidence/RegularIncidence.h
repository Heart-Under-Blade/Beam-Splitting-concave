#pragma once

#include "Incidence.h"

class RegularIncidence : public Incidence
{
public:
	virtual void ComputeDirections(Beam &beam, SplittedBeams<Beam> &beams,
								   Splitting &splitter) const override;
	virtual void ComputeJonesMatrices(Beam &beam, SplittedBeams<Beam> &beams,
									  Splitting &splitter) const override;
};
