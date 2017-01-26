#include "Beam.h"

#include <float.h>
#include <math.h>

#include "macro.h"
#include "geometry_lib.h"
#include <assert.h>

Beam::Beam()
{
	opticalPath = 0;
	e = Point3f(0, 1, 0);
}

void Beam::Copy(const Beam &other)
{
	opticalPath = other.opticalPath;
	D = other.D;
	e = other.e;
	direction = other.direction;

	SetPolygon(other.polygon);

	facetId = other.facetId;
	level = other.level;
	location = other.location;

#ifdef _TRACK_ALLOW
	track = other.track;
#endif
}

Beam::Beam(const Beam &other)
	: JMatrix(other.JMatrix)
{
	Copy(other);
}

Beam::Beam(Beam &&other)
{
	opticalPath = other.opticalPath;
	D = other.D;
	e = other.e;
	direction = other.direction;

	SetPolygon(other.polygon);

	facetId = other.facetId;
	level = other.level;
	location = other.location;

	other.opticalPath = 0;
	other.D = 0;
	other.e = Point3f(0, 0, 0);
	other.direction = Point3f(0, 0, 0);

	other.polygon.size = 0;

	other.facetId = 0;
	other.level = 0;
	other.location = Location::Outside;
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
	const float &y = direction.cy;
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
	if (this != &other)
	{
		Copy(other);
		JMatrix = other.JMatrix;
	}

	return *this;
}

Beam &Beam::operator = (Beam &&other)
{
	if (this != &other)
	{
		opticalPath = other.opticalPath;
		D = other.D;
		e = other.e;
		direction = other.direction;

		SetPolygon(other.polygon);

		facetId = other.facetId;
		level = other.level;
		location = other.location;

		JMatrix = other.JMatrix;

		other.opticalPath = 0;
		other.D = 0;
		other.e = Point3f(0, 0, 0);
		other.direction = Point3f(0, 0, 0);

		other.polygon.size = 0;

		other.facetId = 0;
		other.level = 0;
		other.location = Location::Outside;
	}

	return *this;
}

void Beam::RotatePlane(const Point3f &newBasis)
{
	RotateJMatrix(newBasis);
}

void Beam::RotateJMatrix(const Point3f &newBasis)
{
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

	LOG_ASSERT(fabs(cs) <= 1.0+DBL_EPSILON);

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

void Beam::AddVertex(const Point3f &vertex)
{
	polygon.arr[polygon.size++] = vertex;
}

void Beam::SetPolygon(const Polygon &other)
{
	polygon.size = other.size;

	for (int i = 0; i < other.size; ++i)
	{
		polygon.arr[i] = other.arr[i];
	}
}
