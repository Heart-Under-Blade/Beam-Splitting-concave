#include "Beam.h"

#include <float.h>
#include <math.h>

#include "global.h"
#include "vector_lib.h"

Beam::Beam()
{
	opticalPath = 0;
	shapeSize = 0;
}

void Beam::copy(const Beam &other)
{
	opticalPath = other.opticalPath;
	shapeSize = other.shapeSize;
	D = other.D;
	e = other.e;
	direction = other.direction;

	for (int i = 0; i < other.shapeSize; ++i)
	{
		shape[i] = other.shape[i];
	}
}

Beam::Beam(const Beam &other)
	: JMatrix(other.JMatrix)
{
	copy(other);
}

void Beam::RotateSpherical(const Point3f &dir, const Point3f &polarBasis)
{
	Point3f newBasis;
	double cs = DotProduct(dir, direction);

	if (fabs(1.0 - cs) <= DBL_EPSILON)
	{
		newBasis = -polarBasis;
	}
	else
	{
		if (fabs(1.0 + cs) <= DBL_EPSILON)
		{
			newBasis = polarBasis;
		}
		else
		{
			double fi, teta;
			GetSpherical(fi, teta);
			newBasis = Point3f(-sin(fi),cos(fi), 0);
		}
	}

	RotateJMatrix(newBasis);
}

void Beam::GetSpherical(double &fi, double &teta) const
{
	const float &x = direction.cx;
	const float &y = direction.cx;
	const float &z = direction.cz;

	if (fabs(z + 1.0) < DBL_EPSILON) // forward
	{
		fi = 0;
		teta = M_PI;
		return;
	}

	if (fabs(z - 1.0) < DBL_EPSILON) // bacward
	{
		fi = 0;
		teta = 0;
		return;
	}

	double tmp = y*y;

	if (tmp < DBL_EPSILON)
	{
		tmp = (x > 0) ? 0 : M_PI;
	}
	else
	{
		tmp = acos(x/sqrt(x*x + tmp));

		if (y < 0)
		{
			tmp = M_2PI - tmp;
		}
	}

	fi = (tmp < M_2PI) ? tmp : 0;
	teta = acos(z);
}

Beam & Beam::operator = (const Beam &other)
{
	copy(other);
	JMatrix = other.JMatrix;
	return *this;
}

void Beam::RotatePlane(const Point3f &newBasis)
{
	RotateJMatrix(newBasis);
}

void Beam::RotateJMatrix(const Point3f &newBasis)
{
	/// TODO: потом заменить обратно на DBL_EPSILON
	const double eps = 1e2 * FLT_EPSILON/*DBL_EPSILON*/; // acceptable precision
	double cs = DotProduct(newBasis, e);

	if (fabs(1.0 - cs) < eps)
	{
		return;
	}

	if (fabs(1.0 + cs) < eps)
	{
		JMatrix.m11 = -JMatrix.m11;
		JMatrix.m12 = -JMatrix.m12;
		JMatrix.m21 = -JMatrix.m21;
		JMatrix.m22 = -JMatrix.m22;
		return;
	}

#ifdef _DEBUG
	if(fabs(cs) > 1.0+DBL_EPSILON)  throw " Rot: Error!";
#endif

	Point3f k;
	CrossProduct(e, newBasis, k);
	Normalize(k);

	double angle = acos(cs);

	Point3f r = k + direction;

	if(Norm(r) <= 0.5)
	{
		angle = -angle;
	}

	double sn = sin(angle); // the rotation of matrix "m"

	complex b00 = JMatrix.m11*cs + JMatrix.m21*sn; // first row of the result
	complex b01 = JMatrix.m12*cs + JMatrix.m22*sn;

	JMatrix.m21 = JMatrix.m21*cs - JMatrix.m11*sn;
	JMatrix.m22 = JMatrix.m22*cs - JMatrix.m12*sn;
	JMatrix.m11 = b00;
	JMatrix.m12 = b01;

	e = newBasis;
}

double Beam::Square() const
{
	/// NOTE: не работает на невыпуклых гранях
	double square = 0;

	const Point3f &basePoint = shape[0];
	Point3f p1 = shape[1] - basePoint;

	for (int i = 2; i < shapeSize; ++i)
	{
		Point3f p2 = shape[i] - basePoint;
		Point3f res;
		CrossProduct(p2, p1, res);
		square += sqrt(Norm(res));
		p1 = p2;
	}

	if (square < 0)
	{	/// TODO: для опт. узнать в какую сторону ориентированы точки в пучке
		square *= (-1);
	}

	return square/2.0;
}

Point3f Beam::Center() const
{
	Point3f p(0, 0, 0);

	for (int i = 0; i < shapeSize; ++i)
	{
		p = p + shape[i];
	}

	return p/shapeSize;
}

double Beam::CrossSection() const
{
	const double Eps = 1e7*DBL_EPSILON;

	Point3f normal;
	Point3f p1 = shape[1] - shape[0];
	Point3f p2 = shape[2] - shape[0];
	CrossProduct(p2, p1, normal);

	double e = fabs(DotProduct(normal, direction));

	if (e < Eps)
	{
		return 0;
	}

	double s = Square();
	double n = sqrt(Norm(normal));
	return (e*s) / n; // TODO: опт.
}

void Beam::AddVertex(const Point3f &vertex)
{
	shape[shapeSize] = vertex;
	++shapeSize;
}

void Beam::MulJMatrix(const Beam &other, const complex &coef1, const complex &coef2)
{
	JMatrix.m11 = coef1 * other.JMatrix.m11;
	JMatrix.m12 = coef1 * other.JMatrix.m12;
	JMatrix.m21 = coef2 * other.JMatrix.m21;
	JMatrix.m22 = coef2 * other.JMatrix.m22;
}
