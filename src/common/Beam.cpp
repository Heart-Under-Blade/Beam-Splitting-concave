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

	os << "polygon: {" << endl;
	os << Polygon(beam);
	os << "}" << endl << endl;

	os << "level: " << beam.actNo << endl
	   << "last facet: " << beam.facet->index << endl
	   << "direction: "
	   << beam.direction.coordinates[0] << ", "
	   << beam.direction.coordinates[1] << ", "
	   << beam.direction.coordinates[2] << ", "
	   << beam.direction.d_param << endl << endl;

	return os;
}

Point3d Proj(const Point3d& Tx, const Point3d& Ty, const Point3d& r,  const Point3d& pnt)
{
	const  Point3d p_pr = pnt - r*Point3d::DotProduct(r, pnt); // расчёт коор-т в СК наблюдателя
	return Point3d(Point3d::DotProduct(p_pr, Tx), Point3d::DotProduct(p_pr, Ty), 0); //*/
}

Point3d Proj(const Point3d& _r, const Point3d &pnt)
{
	Point3d _Tx,  // условная горизонталь СК экрана в СК тела
		_Ty;  // третья ось (условная вертикаль СК экрана)
	const double tmp = sqrt(SQR(_r.x)+SQR(_r.y));
	(fabs(_r.z)>1-DBL_EPSILON) ? (_Tx=Point3d(0,-_r.z,0), _Ty=Point3d(1,0,0))
							   : (_Tx=Point3d(_r.y/tmp,-_r.x/tmp,0), _Ty=Point3d::CrossProduct(_r,_Tx));
	return Proj(_Tx, _Ty, _r, pnt);
}

Beam::Beam()
{
	polarizationBasis = Vector3f(0, 1, 0); // OPT: кажется, это используется только для стартового пучка
}

void Beam::Copy(const Beam &other)
{
	Track::operator=(other);

	direction = other.direction;
	polarizationBasis = other.polarizationBasis;
}

Beam::Beam(const Beam &other)
	: Jones(other.Jones)
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

// NOTE: The beam should be always extertal
Vector3f Beam::RotateSpherical(const Vector3f &dir, const Vector3f &polarBasis)
{
	Vector3f newBasis;
	double cs = Point3f::DotProduct(dir, direction);

	if (fabs(1.0 - cs) < DBL_EPSILON)
	{
		newBasis = -polarBasis;
	}
	else
	{
		if (fabs(1.0 + cs) < DBL_EPSILON)
		{
			newBasis = polarBasis;
		}
		else
		{
			Orientation sph = direction.ToOrientation();
			newBasis = Vector3f(-sin(sph.zenith), cos(sph.azimuth), 0);
		}
	}

	RotateJones(newBasis, false);
	return newBasis;
}

Beam &Beam::operator = (const Beam &other)
{
	if (this != &other) // OPT: попробовать убрать это уловие для ускорения
	{
		Copy(other);
		Polygon::operator =(other);
		Jones = other.Jones;
	}

	return *this;
}

Beam &Beam::operator = (const Polygon &other)
{
	Polygon::operator =(other);
	return *this;
}

void Beam::SetDefault(Beam &other)
{
	other.direction = Vector3f(0, 0, 0);
	other.polarizationBasis = Vector3f(0, 0, 0);

	other.facet = nullptr;
	other.actNo = 0;
	other.locations = 0;
	other.id = 0;
}

Beam &Beam::operator = (Beam &&other)
{
	if (this != &other)
	{
		Polygon::operator =(other);
		Copy(other);
		Jones = other.Jones;
		SetDefault(other);
	}

	return *this;
}

void Beam::SetIsInside(bool isIn)
{
	if (!isIn)
	{
		int loc = 1;
		loc <<= actNo;
		locations |= loc;
	}
}

bool Beam::IsShadow() const
{
	return facet->index == INT_MAX;
}

// REF: перенести в Matrix2x2c
void Beam::MultiplyByFresnel(const complex &f1, const complex &f2)
{
	Jones.m11 *= f1;
	Jones.m12 *= f1;
	Jones.m21 *= f2;
	Jones.m22 *= f2;
}

void Beam::RotateJones(const Vector3f &normal, bool isInside)
{
	Point3f newBasis = (isInside) ? Point3f::CrossProduct(normal, direction)
								  : Point3f::CrossProduct(normal, -direction);
	Point3f::Normalize(newBasis);

	const double eps = 1e2 * FLT_EPSILON/*DBL_EPSILON*/; // acceptable precision
	double cs = Point3f::DotProduct(newBasis, polarizationBasis);

	if (fabs(1.0 - cs) >= eps)
	{
		if (fabs(1.0 + cs) < eps)
		{
			Jones.m11 = -Jones.m11;
			Jones.m12 = -Jones.m12;
			Jones.m21 = -Jones.m21;
			Jones.m22 = -Jones.m22;
		}
		else
		{
			Point3f k;
			Point3f::CrossProduct(polarizationBasis, newBasis, k);
			Point3f::Normalize(k);

			double angle = acos(cs);

			Point3f r = k + direction;

			if (Point3f::Norm(r) <= 0.5)
			{
				angle = -angle;
			}

			double sn = sin(angle); // the rotation of matrix "m"

			complex b00 = Jones.m11*cs + Jones.m21*sn; // first row of the result
			complex b01 = Jones.m12*cs + Jones.m22*sn;

			Jones.m21 = Jones.m21*cs - Jones.m11*sn;
			Jones.m22 = Jones.m22*cs - Jones.m12*sn;
			Jones.m11 = b00;
			Jones.m12 = b01;
		}
	}

	polarizationBasis = newBasis;
}

void Beam::SetPolygon(const Polygon &other)
{
	nVertices = other.nVertices;

	for (int i = 0; i < other.nVertices; ++i)
	{
		vertices[i] = other.vertices[i];
	}
}

void Beam::SetDefault()
{
	locations = 0;
	polarizationBasis = Vector3f(0, 1, 0);
}

void Beam::CopyTrack(const Track &other)
{
	Track::operator=(other);
}

