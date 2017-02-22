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

matrixC JonesMatrix::operator *(const complex &z) const
{
	matrixC rez(2, 2);

	rez[0][0] = m11 * z;
	rez[0][1] = m12 * z;
	rez[1][0] = m21 * z;
	rez[1][1] = m22 * z;

	return rez;
}

