#pragma once

#include "Incidence.h"

class RegularIncidence : public Incidence
{
public:
	RegularIncidence(const complex &m_ri);
	void ComputeDirections(Beam &beam, SplittedBeams<Beam> &beams) const;
	void ComputeJonesMatrices(Beam &beam, SplittedBeams<Beam> &beams) const;

	Point3f normal;
};
