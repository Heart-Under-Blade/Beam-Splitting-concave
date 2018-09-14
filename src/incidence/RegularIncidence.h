#pragma once

#include "Incidence.h"

class Beam;
class Splitting;

class RegularIncidence : public Incidence
{
public:
	virtual void ComputeLightParams(const Beam &incidentBeam,
									Splitting &splitter) const override;

	virtual void ComputeJonesMatrices(const Beam &incidentBeam,
									  Splitting &splitter) const override;
};
