#include "JonesMatrix.h"

JonesMatrix::JonesMatrix()
{
	m11 = m22 = 1;
	m12 = m21 = 0;
}

double JonesMatrix::Norm() const
{
	return norm(m11) + norm(m12) + norm(m21) + norm(m22);
}

