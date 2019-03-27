#pragma once

#include "Matrix4x4.h"
#include "JonesMatrix.h"
#include "matrix.hpp"

class MuellerMatrix : public Matrix4x4d
{
public:
	MuellerMatrix();
	MuellerMatrix(const Matrix2x2c &in);
	MuellerMatrix(const matrixC &in);
};
