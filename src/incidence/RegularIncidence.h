#pragma once

#include "Incidence.h"

class RegularIncidence : public Incidence
{
public:
	RegularIncidence();
	void ComputeDirections(Beam &beam, BeamPair<Beam> &beams) const override;
	void ComputeJonesMatrices(Beam &beam, BeamPair<Beam> &beams) const override;
};
