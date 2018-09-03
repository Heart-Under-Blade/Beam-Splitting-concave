#include "Matrix4x4.h"


Matrix4x4d::Matrix4x4d()
{
	Reset();
}

void Matrix4x4d::Reset()
{
	for (int i = 0; i < MX_RANK; ++i)
	{
		for (int j = 0; j < MX_RANK; ++j)
		{
			m_array[i][j] = 0.f;
		}
	}
}

Matrix4x4d Matrix4x4d::operator - (const Matrix4x4d &other) const
{
	Matrix4x4d res;

	for (int i = 0; i < MX_RANK; ++i)
	{
		for (int j = 0; j < MX_RANK; ++j)
		{
			res(i, j) = m_array[i][j] - other(i, j);
		}
	}

	return res;
}

void Matrix4x4d::operator +=(const Matrix4x4d &other)
{
	for (int i = 0; i < MX_RANK; ++i)
	{
		for (int j = 0; j < MX_RANK; ++j)
		{
			m_array[i][j] += other(i, j);
		}
	}
}

void Matrix4x4d::operator *= (double value)
{
	for (int i = 0; i < MX_RANK; ++i)
	{
		for (int j = 0; j < MX_RANK; ++j)
		{
			m_array[i][j] *= value;
		}
	}
}

double Matrix4x4d::operator ()(unsigned int i, unsigned int j) const
{
	return m_array[i][j];
}

double &Matrix4x4d::operator ()(unsigned int i, unsigned int j)
{
	return m_array[i][j];
}

std::ofstream& operator << (std::ofstream& out, const Matrix4x4d& m)
{
	for (int i = 0; i < MX_RANK; i++)
	{
		for (int j = 0; j < MX_RANK; j++)
		{
			out << m(i, j);

			if ((i+1 != MX_RANK) || (j+1 != MX_RANK))
			{
				out << " ";
			}
		}
	}
	return out;
}
