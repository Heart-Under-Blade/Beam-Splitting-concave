#pragma once

#include "matrix.hpp"
#include "compl.hpp"

/**
 * @brief The Matrix2x2c class
 * Squad matrix with 4 complex elements (2x2)
 */
class Matrix2x2c
{
public:
	Matrix2x2c();
	Matrix2x2c(const matrixC &m);

	double Norm() const;
	void Fill(const complex &value);

	complex m11;
	complex m12;
	complex m21;
	complex m22;

	matrixC operator * (const complex &value) const;
	Matrix2x2c &operator *= (const double &value);
	Matrix2x2c &operator += (const Matrix2x2c &other);
};
