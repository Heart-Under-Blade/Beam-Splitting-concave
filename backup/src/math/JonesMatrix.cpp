#include "JonesMatrix.h"

Matrix2x2c::Matrix2x2c()
{
	m11 = m22 = 1;
	m12 = m21 = 0;
}

Matrix2x2c::Matrix2x2c(const matrixC &m)
{
	m11 = m[0][0];
	m12 = m[0][1];
	m21 = m[1][0];
	m22 = m[1][1];
}

double Matrix2x2c::Norm() const
{
	return norm(m11) + norm(m12) + norm(m21) + norm(m22);
}

void Matrix2x2c::Fill(const complex &value)
{
	m12 = m21 = m11 = m22 = value;
}

matrixC Matrix2x2c::operator * (const complex &value) const
{
	matrixC rez(2, 2);

	rez[0][0] = m11 * value;
	rez[0][1] = m12 * value;
	rez[1][0] = m21 * value;
	rez[1][1] = m22 * value;

	return rez;
}

Matrix2x2c &Matrix2x2c::operator *=(const double &value)
{
	m11 *= value;
	m12 *= value;
	m21 *= value;
	m22 *= value;
	return *this;
}

Matrix2x2c &Matrix2x2c::operator +=(const Matrix2x2c &other)
{
	m11 += other.m11;
	m12 += other.m12;
	m21 += other.m21;
	m22 += other.m22;
	return *this;
}

