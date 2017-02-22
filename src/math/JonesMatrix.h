#pragma once

#include "matrix.hpp"
#include "compl.hpp"

class JonesMatrix
{
public:
	JonesMatrix();

	double Norm() const;

	complex m11;
	complex m12;
	complex m21;
	complex m22;

	matrixC operator * (const complex& z) const;
};
