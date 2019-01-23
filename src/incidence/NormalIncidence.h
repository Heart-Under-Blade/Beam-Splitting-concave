#pragma once

#include "RegularIncidence.h"

class NormalIncidence : public Incidence
{
public:
	NormalIncidence();
	void SetSplitting(Splitting *splitting);

	void ComputeDirections(Beam &beam, BeamPair<Beam> &beams);
	void ComputeJonesMatrices(Beam &beam, BeamPair<Beam> &beams);
private:
	complex fresnels[4];
};
