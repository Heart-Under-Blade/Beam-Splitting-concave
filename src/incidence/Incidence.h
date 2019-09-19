#pragma once

#include "Splitting.h"

/**
 * @brief Compute params of splitted beams
 */
class Incidence
{
public:
	Incidence();

	void SetSplitting(Splitting *splitting);

	virtual void ComputeDirections(Beam &beam, BeamPair<Beam> &beams,
								   bool isBeamInside) const = 0;
	virtual void ComputeJonesMatrices(Beam &beam, BeamPair<Beam> &beams,
									  bool isBeamInside) const = 0;

protected:
	Splitting *m_splitting;
};
