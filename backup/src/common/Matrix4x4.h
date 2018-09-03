#pragma once

#include <fstream>

#define MX_RANK 4

class Matrix4x4d
{
public:
	Matrix4x4d();

	void Reset();

	Matrix4x4d operator - (const Matrix4x4d &other) const;
	void operator +=(const Matrix4x4d &other);
	void operator *=(double value);
	double operator ()(unsigned int i, unsigned int j) const;
	double &operator () (unsigned int i, unsigned int j);

	friend std::ofstream& operator << (std::ofstream &out, const Matrix4x4d &m);

protected:
	double m_array[MX_RANK][MX_RANK];
};
