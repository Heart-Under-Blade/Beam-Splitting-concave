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
	   << "location: " << beam.isInside << endl
	   << "direction: "
	   << beam.direction.cx << ", "
	   << beam.direction.cy << ", "
	   << beam.direction.cz << ", "
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
	opticalPath = 0;
	polarizationBasis = Vector3f(0, 1, 0);
}

void Beam::Copy(const Beam &other)
{
	Track::operator=(other);

	opticalPath = other.opticalPath;
	front = other.front;
	direction = other.direction;
	polarizationBasis = other.polarizationBasis;

	isInside = other.isInside;
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
			double phi, teta;
			GetSpherical(phi, teta);
			newBasis = Vector3f(-sin(phi), cos(phi), 0);
		}
	}

	RotateJMatrix(newBasis);
	return newBasis;
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

	other.facet = nullptr;
	other.actNo = 0;
	other.isInside = false;
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
		Jones = other.Jones;
		SetDefault(other);
	}

	return *this;
}

void Beam::SetLocation(bool isIn)
{
	isInside = isIn;

	if (!isInside)
	{	// write location
		int loc = 1;
		loc <<= actNo;
		locations |= loc;
	}
}

// REF: перенести в Matrix2x2c
void Beam::MultiplyJonesMatrix(const complex &f1, const complex &f2)
{
	Jones.m11 *= f1;
	Jones.m12 *= f1;
	Jones.m21 *= f2;
	Jones.m22 *= f2;
}

complex Beam::DiffractionIncline(const Point3d &pt, double wavelength) const
{
	const double eps1 = /*1e9**/100*FLT_EPSILON;
	const double eps2 = /*1e6**/FLT_EPSILON;

	Vector3f _n = Normal();

	int begin, startIndex, endIndex;
	bool order = !(Point3f::DotProduct(_n, direction) < 0);

	if (order)
	{
		begin = 0;
		startIndex = nVertices-1;
		endIndex = -1;
	}
	else
	{
		begin = nVertices-1;
		startIndex = 0;
		endIndex = nVertices;
	}

	Point3d n = Point3d(_n.cx, _n.cy, _n.cz);

	const Point3f &dir = direction;
	Point3d k_k0 = -pt + Point3d(dir.cx, dir.cy, dir.cz);

	Point3f cntr = Center();
	Point3d center = Proj(n, Point3d(cntr.cx, cntr.cy, cntr.cz));

	Point3d	pt_proj = Proj(n, k_k0);

	const double
			A = pt_proj.x,
			B = pt_proj.y;

	complex one(0, -1);
	const double k = M_2PI/wavelength;

	if (fabs(A) < eps2 && fabs(B) < eps2)
	{
		return -one/wavelength*Area();
	}

	complex s(0, 0);

//	std::list<Point3d>::const_iterator p = polygon.arrthis->v.begin();
//	Point3d p1 = Proj(this->N, *p++)-cnt, p2; // переводим вершины в систему координат грани

	Point3d p1 = Proj(n, arr[begin]) - center;
	Point3d p2;

	if (fabs(B) > fabs(A))
	{
		for (int i = startIndex; i != endIndex;)
		{
			p2 = Proj(n, arr[i]) - center;

			if (fabs(p1.x - p2.x) < eps1)
			{
				p1 = p2;

				if (order)
				{
					--i;
				}
				else
				{
					++i;
				}

				continue;
			}

			const double
					ai = (p1.y-p2.y)/(p1.x-p2.x),
					bi = p1.y - ai*p1.x,
					Ci = A+ai*B;

			complex tmp;

			if (fabs(Ci) < eps1)
			{
				tmp = complex(-k*k*Ci*(p2.x*p2.x-p1.x*p1.x)/2.0,k*(p2.x-p1.x));
			}
			else
			{
				tmp = (exp_im(k*Ci*p2.x) - exp_im(k*Ci*p1.x))/Ci; // OPT: (k*Ci)
			}

			s += exp_im(k*B*bi) * tmp;
			p1 = p2;

			if (order)
			{
				--i;
			}
			else
			{
				++i;
			}
		}

		s /= B;
	}
	else
	{
		for (int i = startIndex; i != endIndex;)
		{
			p2 = Proj(n, arr[i]) - center;

			if (fabs(p1.y - p2.y)<eps1)
			{
				p1 = p2;

				if (order)
				{
					--i;
				}
				else
				{
					++i;
				}

				continue;
			}

			const double ci = (p1.x-p2.x)/(p1.y-p2.y),
					di = p1.x - ci*p1.y,
					Ei = A*ci+B;

			s += exp_im(k*A*di) * (fabs(Ei)<eps1 ? complex(-k*k*Ei*(p2.y*p2.y-p1.y*p1.y)/2.0,k*(p2.y-p1.y))
												 : (exp_im(k*Ei*p2.y) - exp_im(k*Ei*p1.y))/Ei);// OPT: (k*Ei)
			p1 = p2;

			if (order)
			{
				--i;
			}
			else
			{
				++i;
			}
		}

		s /= -A;
	}

	return one*wavelength*s/SQR(M_2PI);
}

void Beam::RotateJMatrix(const Vector3f &newBasis)
{
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
}

void Beam::SetPolygon(const Polygon &other)
{
	nVertices = other.nVertices;

	for (int i = 0; i < other.nVertices; ++i)
	{
		arr[i] = other.arr[i];
	}
}

void Beam::Clear()
{
	locations = 0;
	opticalPath = 0;
	polarizationBasis = Vector3f(0, 1, 0);
}

void Beam::AddOpticalPath(double path)
{
	opticalPath += path;
	front = Point3f::DotProduct(-direction, Center());
#ifdef _DEBUG // DEB
	ops.push_back(path);
#endif
}

void Beam::CopyTrack(const Track &other)
{
	Track::operator=(other);
}
