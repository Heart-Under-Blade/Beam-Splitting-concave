#pragma once

#define MAX_FACET_NUM 64

#define MIN_VERTEX_NUM 3
#define MAX_VERTEX_NUM 64

struct IntArray
{
	int arr[64];
	int size = 0;
};

// fast access for Point3f
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

	Point3d(const Point3f &other)
	{
		x = other.cx;
		y = other.cy;
		z = other.cz;
		d = other.d_param;
	}

	Point3d(double p_x, double p_y, double p_z, double p_d = 0.0)
	{
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

struct Polygon
{
	Point3f arr[MAX_VERTEX_NUM];
	int size = 0;
};

struct Facet
{
	Polygon polygon;
	Point3f normal[2]; ///< internal and external
};


/**
 * Functions
 */

float DotProduct(const Point3f &v1, const Point3f &v2);

double DotProductD(const Point3d &v1, const Point3d &v2);

double Norm(const Point3f &point);

void CrossProduct(const Point3f &v1, const Point3f &v2, Point3f &res);

void Normalize(Point3f &v);

Point3f NormalToPolygon(const Point3f *facet);

Point3f CenterOfPolygon(const Polygon &polygon);

double Length(const Point3f &v);
