#pragma once

#include "Incidence.h"

class Beam;
class Splitting;

class NormalIncidence : public Incidence
{
public:
	virtual void ComputeLightParams(const Beam &beam,
									Splitting &splitter) const override;

	virtual void ComputeJonesMatrices(const Beam &beam,
									  Splitting &splitter) const override;
};
