#include "Beam.h"

#include <float.h>
#include <math.h>
#include <assert.h>
#include <list>

#include "macro.h"
#include "geometry_lib.h"

std::ostream& operator << (std::ostream &os, const Beam &beam)
{
	using namespace std;

	os << Polygon(beam);

	os << "level: " << beam.nActs << endl
	   << "last facet: " << beam.lastFacetId << endl
	   << "location: " << beam.location << endl
	   << "direction: "
	   << beam.direction.cx << ", "
	   << beam.direction.cy << ", "
	   << beam.direction.cz << ", "
	   << beam.direction.d_param << endl << endl;

	return os;
}

Beam::Beam()
{
	locations = 0;
	opticalPath = 0;
	polarizationBasis = Vector3f(0, 1, 0);
}

void Beam::Copy(const Beam &other)
{
	opticalPath = other.opticalPath;
	front = other.front;
	direction = other.direction;
	polarizationBasis = other.polarizationBasis;

	lastFacetId = other.lastFacetId;
	nActs = other.nActs;
	location = other.location;
	locations = other.locations;
	id = other.id;
#ifndef _DEBUG // DEB
//	pols = other.pols;
//	ops = other.ops;
//	dirs = other.dirs;
#endif
}

Beam::Beam(const Beam &other)
	: J(other.J)
{
	Copy(other);
	Polygon::operator =(other);
}

Beam::Beam(const Polygon &other)
	: Polygon(other)
{
}

Beam::Beam(Beam &&other)
	: Polygon(other)
{
	Copy(other);
	SetDefault(other);
}

Vector3f Beam::RotateSpherical(const Vector3f &dir, const Vector3f &polarBasis)
{
	Vector3f newBasis;
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
			newBasis = Vector3f(-sin(fi), cos(fi), 0);
		}
	}

	RotateJMatrix(newBasis);
	return newBasis;
}

void Beam::SetMatrix(const Matrix2x2c &matrix)
{
	J = matrix;
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

Beam &Beam::operator = (const Beam &other)
{
	if (this != &other)
	{
		Copy(other);
		Polygon::operator =(other);
		J = other.J;
	}

	return *this;
}

Beam &Beam::operator = (const Polygon &other)
{
	Polygon::operator =(other);
	return *this;
}

Beam &Beam::operator = (const Light &other)
{
	direction = other.direction;
	polarizationBasis = other.polarizationBasis;
	return *this;
}

void Beam::SetDefault(Beam &other)
{
	other.opticalPath = 0;
	other.front = 0;
	other.direction = Vector3f(0, 0, 0);
	other.polarizationBasis = Vector3f(0, 0, 0);

	other.lastFacetId = 0;
	other.nActs = 0;
	other.location = Location::Out;
	other.locations = 0;

#ifdef _TRACK_ALLOW
	other.id = 0;
#endif
}

Beam &Beam::operator = (Beam &&other)
{
	if (this != &other)
	{
		Polygon::operator =(other);
		Copy(other);
		J = other.J;
		SetDefault(other);
	}

	return *this;
}

void Beam::SetTracingParams(int facetId, int actN, Location loc)
{
	lastFacetId = facetId;
	nActs = actN;
	location = loc;

	if (loc == Location::Out)
	{	// write location
		int loc = 1;
		loc <<= nActs;
		locations |= loc;
	}
}

void Beam::MultiplyJonesMatrix(const complex &c1, const complex &c2)
{
	J.m11 *= c1;
	J.m12 *= c1;
	J.m21 *= c2;
	J.m22 *= c2;
}

void Beam::RotateJMatrix(const Vector3f &newBasis)
{
	const double eps = 1e2 * FLT_EPSILON/*DBL_EPSILON*/; // acceptable precision
	double cs = DotProduct(newBasis, polarizationBasis);

	if (fabs(1.0 - cs) >= eps)
	{
		if (fabs(1.0 + cs) < eps)
		{
			J.m11 = -J.m11;
			J.m12 = -J.m12;
			J.m21 = -J.m21;
			J.m22 = -J.m22;
		}
		else
		{
			Point3f k;
			CrossProduct(polarizationBasis, newBasis, k);
			Normalize(k);

			double angle = acos(cs);

			Point3f r = k + direction;

			if (Norm(r) <= 0.5)
			{
				angle = -angle;
			}

			double sn = sin(angle); // the rotation of matrix "m"

			complex b00 = J.m11*cs + J.m21*sn; // first row of the result
			complex b01 = J.m12*cs + J.m22*sn;

			J.m21 = J.m21*cs - J.m11*sn;
			J.m22 = J.m22*cs - J.m12*sn;
			J.m11 = b00;
			J.m12 = b01;
		}
	}
}

void Beam::SetPolygon(const Polygon &other)
{
	nVertices = other.nVertices;

	for (int i = 0; i < other.nVertices; ++i)
	{
		arr[i] = other.arr[i];
	}
}

void Beam::SetLight(const Vector3f &dir, const Vector3f &polarBasis)
{
	direction = dir;
	polarizationBasis = polarBasis;
}

void Beam::SetLight(const Light &other)
{
	direction = other.direction;
	polarizationBasis = other.polarizationBasis;
}

void Beam::AddOpticalPath(double path)
{
	opticalPath += path;
	front = DotProduct(-direction, Center());
}

void Beam::CopyTrack(const Track &other)
{
	id = other.id;
	locations = other.locations;
}
