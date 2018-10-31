#pragma once

#include "Beam.h"
#include "Splitting.h"

class Incidence
{
public:
	virtual void ComputeDirections(Beam &/*beam*/, Splitting &/*splitter*/) const {};
	virtual void ComputeJonesMatrices(Beam &/*beam*/, Splitting &/*splitter*/) const {};

    void ComputeOpticalPaths(const Beam &beam, Splitting &splitter) const;
};
