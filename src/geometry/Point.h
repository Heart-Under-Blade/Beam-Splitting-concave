#pragma once

#include <math.h>
#include <iostream>

// short access for Point3f
#define d_param coordinates[3]

/**
 * @brief A 3-dimentional point
 */
struct Point3f
{
	float coordinates[4];

    Point3f() {}
    Point3f(float x, float y, float z);
    Point3f(float x, float y, float z, float d);
    Point3f(const Point3f &other);

	bool IsEqualTo(const Point3f &other, float eps) const;

    Point3f &operator = (const Point3f &other);
    Point3f operator * (double value) const;
	Point3f operator *= (double value);
    Point3f operator / (double value) const;
    Point3f operator - (const Point3f &value) const;
    Point3f operator + (const Point3f &value) const;
    Point3f operator += (double value);
	Point3f operator - () const;

	friend std::ostream & operator << (std::ostream &os, const Point3f &p);

    static float DotProduct(const Point3f &v1, const Point3f &v2);
    static void CrossProduct(const Point3f &v1, const Point3f &v2, Point3f &res);
    static Point3f CrossProduct(const Point3f &v1, const Point3f &v2);
	static double Norm(const Point3f &v);
    static void Normalize(Point3f &v);
	static double Length(const Point3f &v);

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
		x = other.coordinates[0];
		y = other.coordinates[1];
		z = other.coordinates[2];
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


    static double DotProduct(const Point3d &v1, const Point3d &v2)
    {
        return	  v1.x * v2.x
                + v1.y * v2.y
                + v1.z * v2.z;
    }

    static double Norm(const Point3d &p)
    {
        return	  p.x * p.x
                + p.y * p.y
                + p.z * p.z;
    }

    static Point3d CrossProduct(const Point3d &v1, const Point3d &v2)
    {
        Point3d res;
        res.x = v1.y*v2.z - v1.z*v2.y;
        res.y = v1.z*v2.x - v1.x*v2.z;
        res.z = v1.x*v2.y - v1.y*v2.x;
        return res;
    }

    static double Length(const Point3d &v)
    {
        return sqrt(Norm(v));
    }

    static Point3d Normalize(const Point3d &v)
    {
        Point3d res;
        double lenght = sqrt(Norm(v));
        res.x = v.x / lenght;
        res.y = v.y / lenght;
        res.z = v.z / lenght;
        return res;
    }
};
