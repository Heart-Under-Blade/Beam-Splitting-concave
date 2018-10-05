#pragma once

#include "Incidence.h"

class Beam;
class Splitting;

class RegularIncidence : public Incidence
{
public:
	virtual void ComputeDirections(const Beam &beam,
									Splitting &splitter) const override;

	virtual void ComputeJonesMatrices(const Beam &beam,
									  Splitting &splitter) const override;
};
