#pragma once

#include "Beam.h"

class PhisBeam : public Beam
{
public:
	PhisBeam();

public:
	double opticalPath;				///< optical path of beam
	double frontPosition;						///< current position of phase front from Ax+By+Cz+D=0
};
