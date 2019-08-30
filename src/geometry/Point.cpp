#include "Point.h"
#include "intrinsic/intrinsics.h"
#include "float.h"
#include "common.h"

Point3f::Point3f(float x, float y, float z)
{
	coordinates[0] = x;
	coordinates[1] = y;
	coordinates[2] = z;
}

Point3f::Point3f(float x, float y, float z, float d)
{
	coordinates[0] = x;
	coordinates[1] = y;
	coordinates[2] = z;
	coordinates[3] = d;
}

Point3f::Point3f(const Point3f &other)
{
	coordinates[0] = other.coordinates[0];
	coordinates[1] = other.coordinates[1];
	coordinates[2] = other.coordinates[2];
}

Point3f &Point3f::operator = (const Point3f &other)
{
	coordinates[0] = other.coordinates[0];
	coordinates[1] = other.coordinates[1];
	coordinates[2] = other.coordinates[2];

	return *this;
}

Point3f Point3f::operator * (double value) const
{
	return Point3f(coordinates[0] * value,
			coordinates[1] * value,
			coordinates[2] * value);
}

Point3f Point3f::operator *=(double value)
{
	*this = operator*(value);
	return *this;
}

Point3f Point3f::operator / (double value) const
{
	return Point3f(coordinates[0] / value,
			coordinates[1] / value,
			coordinates[2] / value);
}

Point3f Point3f::operator /=(double value)
{
	*this = operator/(value);
	return *this;
}

Point3f Point3f::operator + (const Point3f &value) const
{
	return Point3f(coordinates[0] + value.coordinates[0],
			coordinates[1] + value.coordinates[1],
			coordinates[2] + value.coordinates[2]);
}

Point3f Point3f::operator +=(const Point3f &value)
{
	*this = operator+(value);
	return *this;
}

Point3f Point3f::operator += (double value)
{
	return *this = Point3f(coordinates[0] + value,
			coordinates[1] + value,
			coordinates[2] + value);
}

Point3f Point3f::operator - () const
{
	return Point3f(-coordinates[0], -coordinates[1], -coordinates[2],
			-coordinates[3]);
}

Orientation Point3f::ToOrientation() const
{
	Orientation orient;

	const float &x = coordinates[0];
	const float &y = coordinates[1];
	const float &z = coordinates[2];

	if (fabs(z + 1.0) < DBL_EPSILON) // forward
	{
		orient.zenith = 0;
		orient.azimuth = M_PI;
		return orient;
	}

	if (fabs(z - 1.0) < DBL_EPSILON) // backward
	{
		orient.zenith = 0;
		orient.azimuth = 0;
		return orient;
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

	orient.zenith = (tmp < M_2PI) ? tmp : 0;
	orient.azimuth = acos(z);

	return orient;
}

std::ostream &operator <<(std::ostream &os, const Point3f &p)
{
	os << p.coordinates[0] << " " << p.coordinates[1] << " " << p.coordinates[2];
	return os;
}

float Point3f::DotProduct(const Point3f &v1, const Point3f &v2)
{
	__m128 _v1 = _mm_setr_ps(
				v1.coordinates[0], v1.coordinates[1], v1.coordinates[2], 0.0);
	__m128 _v2 = _mm_setr_ps(
				v2.coordinates[0], v2.coordinates[1], v2.coordinates[2], 0.0);
	__m128 _dp0 = _mm_dp_ps(_v1, _v2, MASK_FULL);
	return _dp0[0];
}

void Point3f::CrossProduct(const Point3f &v1, const Point3f &v2, Point3f &res)
{
	__m128 _v1 = _mm_setr_ps(
				v1.coordinates[0], v1.coordinates[1], v1.coordinates[2], 0.0);
	__m128 _v2 = _mm_setr_ps(
				v2.coordinates[0], v2.coordinates[1], v2.coordinates[2], 0.0);
	__m128 _cp = _cross_product(_v1, _v2);

	res.coordinates[0] = _cp[0];
	res.coordinates[1] = _cp[1];
	res.coordinates[2] = _cp[2];
}

Point3f Point3f::CrossProduct(const Point3f &v1, const Point3f &v2)
{
	__m128 _v1 = _mm_setr_ps(
				v1.coordinates[0], v1.coordinates[1], v1.coordinates[2], 0.0);
	__m128 _v2 = _mm_setr_ps(
				v2.coordinates[0], v2.coordinates[1], v2.coordinates[2], 0.0);
	__m128 _cp = _cross_product(_v1, _v2);

	Point3f res;
	res.coordinates[0] = _cp[0];
	res.coordinates[1] = _cp[1];
	res.coordinates[2] = _cp[2];

	return res;
}

double Point3f::Norm(const Point3f &v)
{
	return	  v.coordinates[0] * v.coordinates[0]
			+ v.coordinates[1] * v.coordinates[1]
			+ v.coordinates[2] * v.coordinates[2];
}

void Point3f::Normalize(Point3f &v)
{
	double lenght = Length(v);
	v.coordinates[0] /= lenght;
	v.coordinates[1] /= lenght;
	v.coordinates[2] /= lenght;
}

double Point3f::Length(const Point3f &v)
{
	return sqrt(Norm(v));
}

bool Point3f::IsEqualTo(const Point3f &other, float eps) const
{
	return (fabs(coordinates[0] - other.coordinates[0]) +
			fabs(coordinates[1] - other.coordinates[1]) +
			fabs(coordinates[2] - other.coordinates[2]))/3 < eps;
}

double Point3f::AbsSum() const
{
	return fabs(coordinates[0]) + fabs(coordinates[1]) + fabs(coordinates[2]);
}

Point3f Point3f::operator - (const Point3f &value) const
{
	return Point3f(coordinates[0] - value.coordinates[0],
			coordinates[1] - value.coordinates[1],
			coordinates[2] - value.coordinates[2]);
}
