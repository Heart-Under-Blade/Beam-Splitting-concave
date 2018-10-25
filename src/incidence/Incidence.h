#pragma once

#include "Beam.h"

class Splitting;

class Incidence
{
public:
	virtual void ComputeDirections(Beam &/*beam*/, Splitting &/*splitter*/) const {};
	virtual void ComputeJonesMatrices(Beam &/*beam*/, Splitting &/*splitter*/) const {};

    void ComputeOpticalPaths(const Beam &beam, Splitting &splitter) const;
};
