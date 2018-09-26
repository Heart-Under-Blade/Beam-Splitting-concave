#include "Point.h"
#include "intrinsic/intrinsics.h"

Point3f::Point3f(float x, float y, float z)
{
    point[0] = x;
    point[1] = y;
    point[2] = z;
}

Point3f::Point3f(float x, float y, float z, float d)
{
    point[0] = x;
    point[1] = y;
    point[2] = z;
    point[3] = d;
}

Point3f::Point3f(const Point3f &other)
{
    point[0] = other.point[0];
    point[1] = other.point[1];
    point[2] = other.point[2];
}

Point3f &Point3f::operator = (const Point3f &other)
{
    point[0] = other.point[0];
    point[1] = other.point[1];
    point[2] = other.point[2];

    return *this;
}

Point3f Point3f::operator * (double value) const
{
    return Point3f(point[0] * value,
            point[1] * value,
            point[2] * value);
}

Point3f Point3f::operator / (double value) const
{
    return Point3f(point[0] / value,
            point[1] / value,
            point[2] / value);
}

Point3f Point3f::operator + (const Point3f &value) const
{
    return Point3f(point[0] + value.point[0],
            point[1] + value.point[1],
            point[2] + value.point[2]);
}

Point3f Point3f::operator += (double value)
{
    return *this = Point3f(point[0] + value,
            point[1] + value,
            point[2] + value);
}

Point3f Point3f::operator - () const
{
	return Point3f(-point[0], -point[1], -point[2], -point[3]);
}

float Point3f::DotProduct(const Point3f &v1, const Point3f &v2)
{
	__m128 _v1 = _mm_setr_ps(v1.cx, v1.cy, v1.cz, 0.0);
	__m128 _v2 = _mm_setr_ps(v2.cx, v2.cy, v2.cz, 0.0);
	__m128 _dp0 = _mm_dp_ps(_v1, _v2, MASK_FULL);
	return _dp0[0];
}

void Point3f::CrossProduct(const Point3f &v1, const Point3f &v2, Point3f &res)
{
	__m128 _v1 = _mm_setr_ps(v1.cx, v1.cy, v1.cz, 0.0);
	__m128 _v2 = _mm_setr_ps(v2.cx, v2.cy, v2.cz, 0.0);
	__m128 _cp = _cross_product(_v1, _v2);

	res.cx = _cp[0];
	res.cy = _cp[1];
	res.cz = _cp[2];
}

Point3f Point3f::CrossProduct(const Point3f &v1, const Point3f &v2)
{
	__m128 _v1 = _mm_setr_ps(v1.cx, v1.cy, v1.cz, 0.0);
	__m128 _v2 = _mm_setr_ps(v2.cx, v2.cy, v2.cz, 0.0);
	__m128 _cp = _cross_product(_v1, _v2);

	Point3f res;
	res.cx = _cp[0];
	res.cy = _cp[1];
	res.cz = _cp[2];

	return res;
}

double Point3f::Norm(const Point3f &point)
{
	return	  point.cx * point.cx
			+ point.cy * point.cy
			+ point.cz * point.cz;
}

void Point3f::Normalize(Point3f &v)
{
	double lenght = sqrt(Norm(v));
	v.cx /= lenght;
	v.cy /= lenght;
	v.cz /= lenght;
}

double Point3f::Length(const Point3f &v)
{
	return sqrt(Norm(v));
}

Point3f Point3f::operator - (const Point3f &value) const
{
    return Point3f(point[0] - value.point[0],
            point[1] - value.point[1],
            point[2] - value.point[2]);
}
