#pragma once

#include "PathedBeam.h"
#include "Splitting.h"

/**
 * @brief Compute params of splitted beams
 */
class Incidence
{
public:
	Incidence();

	void SetSplitting(Splitting *splitting);

	virtual void ComputeDirections(Beam &/*beam*/, BeamPair<Beam> &/*beams*/) const {};
	virtual void ComputeJonesMatrices(Beam &/*beam*/, BeamPair<Beam> &/*beams*/) const {};

	virtual void ComputeOpticalPaths(const PathedBeam &beam,
									 BeamPair<PathedBeam> &beams) const;

protected:
	Splitting *m_splitting;
};
