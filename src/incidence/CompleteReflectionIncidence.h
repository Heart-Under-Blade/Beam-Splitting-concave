#pragma once

#include "Incidence.h"

class Beam;
class Splitting;

class CompleteReflectionIncidence : public Incidence
{
public:
	virtual void ComputeDirections(const Beam &incidentBeam,
									Splitting &splitter) const override;

	virtual void ComputeJonesMatrices(const Beam &incidentBeam,
									  Splitting &splitter) const override;

	void ComputeOpticalPaths(const Beam &incidentBeam,
							 Splitting &splitter) const;
};
