#pragma once

#include "RegularIncidence.h"

class NormalIncidence : public Incidence
{
public:
	NormalIncidence(const complex &ri);

	void ComputeDirections(Beam &beam, SplittedBeams<Beam> &beams);
	void ComputeJonesMatrices(Beam &beam, SplittedBeams<Beam> &beams);
private:
	complex fresnels[4];
};
