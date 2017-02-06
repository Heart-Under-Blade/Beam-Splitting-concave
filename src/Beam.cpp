#include "Beam.h"

#include <float.h>
#include <math.h>
#include <assert.h>
#include <list>

#include "macro.h"
#include "geometry_lib.h"

Point3d Proj(const Point3d& Tx, const Point3d& Ty, const Point3d& r,  const Point3d& pnt)
{
	const  Point3d p_pr = pnt - r*DotProductD(r, pnt); // расчёт коор-т в СК наблюдателя
	return Point3d(DotProductD(p_pr, Tx), DotProductD(p_pr, Ty), 0); //*/
}

Point3d Proj(const Point3d& _r, const Point3d &pnt)
{
	Point3d _Tx,  // условная горизонталь СК экрана в СК тела
		_Ty;  // третья ось (условная вертикаль СК экрана)
	const double tmp = sqrt(SQR(_r.x)+SQR(_r.y));
	(fabs(_r.z)>1-DBL_EPSILON) ? (_Tx=Point3d(0,-_r.z,0), _Ty=Point3d(1,0,0))
							   : (_Tx=Point3d(_r.y/tmp,-_r.x/tmp,0), _Ty=CrossProductD(_r,_Tx));
	return Proj(_Tx, _Ty, _r, pnt);
}

//Point2D Proj(const Point3D& a, int i){
//		return (i==0) ? Point2D(a.y, a.z) :
//		( (i==1) ? Point2D(a.x, a.z) :
//		Point2D(a.x, a.y) );
//	}

Beam::Beam()
{
	opticalPath = 0;
	e = Point3f(0, 1, 0);
}

void Beam::Copy(const Beam &other)
{
//	T = other.T;
//	F = other.F;
//	N = other.N;
	opticalPath = other.opticalPath;
	D = other.D;
	e = other.e;
	direction = other.direction;

	SetPolygon(other.polygon);

	facetID = other.facetID;
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
//	T = other.T;
//	F = other.F;
//	N = other.N;
	opticalPath = other.opticalPath;
	D = other.D;
	e = other.e;
	direction = other.direction;

	SetPolygon(other.polygon);

	facetID = other.facetID;
	level = other.level;
	location = other.location;

	other.opticalPath = 0;
	other.D = 0;
	other.e = Point3f(0, 0, 0);
	other.direction = Point3f(0, 0, 0);

	other.polygon.size = 0;

	other.facetID = 0;
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
//		T = other.T;
//		F = other.F;
//		N = other.N;
		opticalPath = other.opticalPath;
		D = other.D;
		e = other.e;
		direction = other.direction;

		SetPolygon(other.polygon);

		facetID = other.facetID;
		level = other.level;
		location = other.location;

		JMatrix = other.JMatrix;

		other.opticalPath = 0;
		other.D = 0;
		other.e = Point3f(0, 0, 0);
		other.direction = Point3f(0, 0, 0);

		other.polygon.size = 0;

		other.facetID = 0;
		other.level = 0;
		other.location = Location::Outside;
	}

	return *this;
}

complex Beam::DiffractionIncline(const Point3d &pt, double lam) const
{
	const double eps1 = 1e9*DBL_EPSILON;
	const double eps2 = 1e6*DBL_EPSILON;

	Point3d k_k0 = -pt + direction;
	Point3f center = CenterOfPolygon(polygon);

	Point3f n = NormalToPolygon(polygon);
	Point3d	pt_proj = Proj(Point3d(n.cx, n.cy, n.cz), k_k0);
//	Point3d	center = Proj(this->N, r0);

	const double
			A = pt_proj.x,
			B = pt_proj.y;

	complex one(0, -1);
	const double k = M_2PI/lam;

	if (fabs(A) < eps2 && fabs(B) < eps2)
	{
		return -one/lam*AreaOfPolygon(polygon);
	}

	complex s(0, 0);

//	std::list<Point3d>::const_iterator p = polygon.arrthis->v.begin();
//	Point3d p1 = Proj(this->N, *p++)-cnt, p2; // переводим вершины в систему координат грани

	Point3d p1 = polygon.arr[polygon.size-1] - center;
	Point3d p2;

	if (fabs(B) > fabs(A))
	{
		for (unsigned int i = 0; i < polygon.size; ++i)
		{
			p2 = polygon.arr[i] - center;

			if (fabs(p1.x - p2.x) < eps1)
			{
				p1 = p2;
				continue;
			}

			const double
					ai = (p1.y-p2.y)/(p1.x-p2.x),
					bi = p1.y - ai*p1.x,
					Ci = A+ai*B;

			s += exp_im(k*B*bi)* (fabs(Ci) < eps1 ? complex(-k*k*Ci*(p2.x*p2.x-p1.x*p1.x)/2.0,k*(p2.x-p1.x))
												  : (exp_im(k*Ci*p2.x) - exp_im(k*Ci*p1.x))/Ci);
			p1 = p2;
		}

		s /= B;
	}
	else
	{
		for (unsigned int i = 0; i < polygon.size; ++i)
		{
			p2 = polygon.arr[i] - center;

			if (fabs(p1.y - p2.y)<eps1)
			{
				p1 = p2; continue;
			}

			const double ci = (p1.x-p2.x)/(p1.y-p2.y),
					di = p1.x - ci*p1.y,
					Ei = A*ci+B;

			s += exp_im(k*A*di) * (fabs(Ei)<eps1 ? complex(-k*k*Ei*(p2.y*p2.y-p1.y*p1.y)/2.0,k*(p2.y-p1.y))
												 : (exp_im(k*Ei*p2.y) - exp_im(k*Ei*p1.y))/Ei);
			p1 = p2;
		}

		s /= -A;
	}

	return one*lam*s/SQR(M_2PI);
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
