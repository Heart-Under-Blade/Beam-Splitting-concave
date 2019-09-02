#pragma once

#include <ostream>
#include <math.h>

#define CLIP_RESULT_SINGLE 1

#define MIN_VERTEX_NUM 3		///< minimum number of vertices in polygon
#define MAX_VERTEX_NUM 64		///< maximum number of vertices in polygon
#define MAX_POLYGON_NUM 512		///< maximum number of polygons in array of polygons
#define MAX_FACET_NUM 256

class IntArray
{
public:
	int arr[MAX_FACET_NUM];
	size_t size = 0;

	void Add(int elem)
	{
		arr[size++] = elem;
	}
};

template <class T>
class Couple
{
public:
	T first;
	T last;
};

// short access for Point3f
#define cx		coordinates[0]
#define cy		coordinates[1]
#define cz		coordinates[2]
#define d_param coordinates[3]

// short access for normals of Facet
#define in_normal normal[0]
#define ex_normal normal[1]

template <class T>
class Array
{
public:
	T elems[MAX_FACET_NUM];
	int nElems = 0;

	void Add(T elem)
	{
		elems[nElems++] = elem;
	}
};

struct Point2f
{
	float x;
	float y;

	Point2f() {}

	Point2f(float x, float y)
		: x(x), y(y)
	{
	}

	Point2f(const Point2f &other)
	{
		x = other.x;
		y = other.y;
	}

	bool IsEqualTo(const Point2f &other, float eps) const
	{
		return fabs(x - other.x) < eps &&
				fabs(y - other.y) < eps;
	}

	Point2f operator - (const Point2f &other) const
	{
		return Point2f(x - other.x,
					   y - other.y);
	}

	Point2f &operator = (const Point2f &other)
	{
		x = other.x;
		y = other.y;
		return *this;
	}

	double CrossProduct(const Point2f &other)
	{
		return x*other.y - y*other.x;
	}
};

/**
 * @brief The Point3 struct
 * 3-coordinated point
 */
struct Point3f
{
	float coordinates[4]; /// coordinates

	Point3f() {}

	Point3f(float x, float y, float z)
	{
		coordinates[0] = x;
		coordinates[1] = y;
		coordinates[2] = z;
	}

	Point3f(float x, float y, float z, float d)
	{
		coordinates[0] = x;
		coordinates[1] = y;
		coordinates[2] = z;
		coordinates[3] = d;
	}

	Point3f(const Point2f &other)
	{
		coordinates[0] = other.x;
		coordinates[1] = other.y;
		coordinates[2] = 0;
	}

	bool IsEqualTo(const Point3f &other, float eps) const
	{
		return (fabs(coordinates[0] - other.coordinates[0]) +
				fabs(coordinates[1] - other.coordinates[1]) +
				fabs(coordinates[2] - other.coordinates[2]))/3 < eps;
	}

	void Clone(const Point3f &other)
	{
		coordinates[0] = other.coordinates[0];
		coordinates[1] = other.coordinates[1];
		coordinates[2] = other.coordinates[2];
		coordinates[3] = other.coordinates[3];
	}

	friend std::ostream & operator << (std::ostream &os, const Point3f &p);

	Point3f(const Point3f &other)
	{
		coordinates[0] = other.coordinates[0];
		coordinates[1] = other.coordinates[1];
		coordinates[2] = other.coordinates[2];
	}

	Point3f & operator = (const Point3f &other)
	{
		coordinates[0] = other.coordinates[0];
		coordinates[1] = other.coordinates[1];
		coordinates[2] = other.coordinates[2];

		return *this;
	}

	Point3f operator * (double value) const
	{
		return Point3f(coordinates[0] * value,
				coordinates[1] * value,
				coordinates[2] * value);
	}

	Point3f operator / (double value) const
	{
		return Point3f(coordinates[0] / value,
				coordinates[1] / value,
				coordinates[2] / value);
	}

	Point3f operator - (const Point3f &value) const
	{
		return Point3f(coordinates[0] - value.coordinates[0],
				coordinates[1] - value.coordinates[1],
				coordinates[2] - value.coordinates[2]);
	}

	Point3f operator + (const Point3f &value) const
	{
		return Point3f(coordinates[0] + value.coordinates[0],
				coordinates[1] + value.coordinates[1],
				coordinates[2] + value.coordinates[2]);
	}

	Point3f operator += (double value)
	{
		return *this = Point3f(coordinates[0] + value,
				coordinates[1] + value,
				coordinates[2] + value);
	}

	Point3f operator - () const
	{
		return Point3f(-coordinates[0], -coordinates[1], -coordinates[2], -coordinates[3]);
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

	Point3d operator + (const Point3d &other) const
	{
		return Point3d(x + other.x, y + other.y, z + other.z);
	}

	Point3d operator - (const Point3d &other) const
	{
		return Point3d(x - other.x, y - other.y, z - other.z);
	}

	Point3d operator - () const
	{
		return Point3d(-x, -y, -z);
	}

	Point3d operator / (double value) const
	{
		return Point3d(x/value, y/value, z/value);
	}
};

typedef Point3f Vector3f;
typedef Point3d Vector3d;

/**
 * Functions
 */

float DotProduct(const Vector3f &v1, const Vector3f &v2);
double DotProductD(const Vector3d &v1, const Vector3d &v2);
void CrossProduct(const Vector3f &v1, const Vector3f &v2, Vector3f &res);
Point3f CrossProduct(const Point3f &v1, const Point3f &v2);
Point3f IntersectVectors(const Point3f &c1, const Point3f &c2,
						 const Point3f &v1, const Point3f &v2,
						 const Point3f &normalToFacet, bool &isOk);

double Norm(const Vector3f &point);
void Normalize(Vector3f &v);
Vector3d NormalizeD(const Vector3d &v);
double Length(const Vector3f &v);

Point3f ProjectPointToPlane(const Point3f &point, const Vector3f &direction,
							const Vector3f &planeNormal);

// REF: try to create template class 'Array<type>'

/**
 * Functions
 */

double NormD(const Vector3d &point);
Point3d CrossProductD(const Point3d &v1, const Point3d &v2);
double LengthD(const Vector3d &v);
