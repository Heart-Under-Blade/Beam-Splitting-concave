#include "MullerMatrix.h"
#include "compl.hpp"

MuellerMatrix::MuellerMatrix()
	: Matrix4x4d()
{

}

MuellerMatrix::MuellerMatrix(const Matrix2x2c &in)
{
	const double a11 = norm(in.m11), a12 = norm(in.m12),
				 a21 = norm(in.m21), a22 = norm(in.m22);

	double A1 = a11+a21, A2 = a12+a22;

	m_array[0][0] = (A1+A2)/2.0;
	m_array[0][1] = (A1-A2)/2.0;

	A1 = a11-a21;
	A2 = a12-a22;

	m_array[1][0] = (A1+A2)/2.0;
	m_array[1][1] = (A1-A2)/2.0;

	complex C1 = in.m11*conj(in.m12), C2 = in.m22*conj(in.m21);

	m_array[0][2] = -real(C1)-real(C2);
	m_array[0][3] = imag(C2)-imag(C1);
	m_array[1][2] = real(C2)-real(C1);
	m_array[1][3] = -imag(C1)-imag(C2);

	C1 = in.m11*conj(in.m21);
	C2 = in.m22*conj(in.m12);

	m_array[2][0] = -real(C1)-real(C2);
	m_array[2][1] = real(C2)-real(C1);
	m_array[3][0] = imag(C1)-imag(C2);
	m_array[3][1] = imag(C2)+imag(C1);

	C1 = in.m11*conj(in.m22);
	C2 = in.m12*conj(in.m21);

	m_array[2][2] = real(C1)+real(C2);
	m_array[2][3] = imag(C1)-imag(C2);
	m_array[3][2] = -imag(C1)-imag(C2);
	m_array[3][3] = real(C1)-real(C2);
}
