#pragma once

#include "Incidence.h"

class TotalReflectionIncidence : public Incidence
{
public:
	TotalReflectionIncidence();
	void ComputeDirections(Beam &beam, BeamPair<Beam> &beams,
						   bool isBeamInside) const override;
	void ComputeJonesMatrices(Beam &parentBeam, BeamPair<Beam> &beams,
							  bool isBeamInside) const override;
};
