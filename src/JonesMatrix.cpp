#include "JonesMatrix.h"

JonesMatrix::JonesMatrix()
{
	m11 = 0;
	m22 = 1;
}

double JonesMatrix::Norm() const
{
	return norm(m11) + norm(m12) + norm(m21) + norm(m22);
}

