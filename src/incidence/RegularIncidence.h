#pragma once

#include "Incidence.h"

class Beam;
class Splitting;

class RegularIncidence : public Incidence
{
public:
	virtual void ComputeDirections(Beam &beam, Splitting &splitter) const override;
	virtual void ComputeJonesMatrices(Beam &beam, Splitting &splitter) const override;
};
