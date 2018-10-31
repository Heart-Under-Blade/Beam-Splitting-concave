#pragma once

#include "RegularIncidence.h"

class NormalIncidence : public RegularIncidence
{
public:
	static void ComputeDirections(Beam &beam, Splitting &splitter);
	static void ComputeJonesMatrices(Beam &beam, Splitting &splitter);
};
