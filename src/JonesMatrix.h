#pragma once

#include "math/compl.hpp"

class JonesMatrix
{
public:
	JonesMatrix();

	double Norm() const;

	complex m11;
	complex m12;
	complex m21;
	complex m22;
};
