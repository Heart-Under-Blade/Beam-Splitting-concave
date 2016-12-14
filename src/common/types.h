#pragma once

/// fast access for Point3f
#define cx point[0]
#define cy point[1]
#define cz point[2]
#define d_param point[3]

/**
 * @brief The Point3 struct
 * 3D coordinate point
 */
struct Point3f
{
	float point[4]; /// coordinates

	Point3f() {}

	Point3f(float x, float y, float z)
	{
		point[0] = x;
		point[1] = y;
		point[2] = z;
	}

	Point3f(const Point3f &other)
	{
		point[0] = other.point[0];
		point[1] = other.point[1];
		point[2] = other.point[2];
	}

	Point3f & operator = (const Point3f &other)
	{
		point[0] = other.point[0];
		point[1] = other.point[1];
		point[2] = other.point[2];

		return *this;
	}

	Point3f operator * (double value) const
	{
		return Point3f(point[0] * value,
				point[1] * value,
				point[2] * value);
	}

	Point3f operator / (double value) const
	{
		return Point3f(point[0] / value,
				point[1] / value,
				point[2] / value);
	}

	Point3f operator - (const Point3f &value) const
	{
		return Point3f(point[0] - value.point[0],
				point[1] - value.point[1],
				point[2] - value.point[2]);
	}

	Point3f operator + (const Point3f &value) const
	{
		return Point3f(point[0] + value.point[0],
				point[1] + value.point[1],
				point[2] + value.point[2]);
	}

	Point3f operator - () const
	{
		return Point3f(-point[0], -point[1], -point[2]);
	}

} __attribute__ ((aligned (16)));

struct Point3d
{
	double x;
	double y;
	double z;
	double d;

	Point3d() {}

	Point3d(double p_x, double p_y, double p_z, double p_d = 0.0) {
		x = p_x;
		y = p_y;
		z = p_z;
		d = p_d;
	}

	Point3d operator * (double value) const
	{
		return Point3d(x*value, y*value, z*value);
	}

	Point3d operator - (const Point3d &value) const
	{
		return Point3d(x - value.x, y - value.y, z - value.z);
	}
};

/**
 * @brief The Facet struct
 * Facet of polyhedron
 */
//struct Facet
//{
//	Point3f vertices[32];
//	Point3f normal;
//	int size;

//	Facet() {}
//	Facet(int p_size) {
//		size = p_size;
//	}

//} __attribute__ ((aligned (16)));

struct Polygon
{
	Point3f points[32];
	int size;

	Polygon() {}
	Polygon(int p_size) {
		size = p_size;
	}

} __attribute__ ((aligned (16)));
