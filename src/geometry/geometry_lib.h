#pragma once

#include "Polygon.h"
#include "Intersection.h"
#include <ostream>

#define CLIP_RESULT_SINGLE 1

#define MAX_FACET_NUM 256
#define ROT_MTR_RANK 3

class Facet;

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

template <class T>
class Couple
{
public:
	T first;
	T last;
};

/**
 * The units (degrees or radians)
 * is defined by yours i.e. if you consider units of this object
 * as degrees from beginning then you can call ToRadian() to convert angles
 * of the object from degrees to radians.
 */
class Angle3d
{
public:
	double alpha;
	double beta;
	double gamma;

	Angle3d() {}
	Angle3d(double a, double b, double g)
	{
		alpha = a;
		beta = b;
		gamma = g;
	}

	/**
	 * @brief Convert angles to radians
	 */
	void ToRadian()
	{
		alpha = DegToRad(alpha);
		beta = DegToRad(beta);
		gamma = DegToRad(gamma);
	}

	Angle3d ToRadian() const
	{
		Angle3d angle = *this;
		angle.ToRadian();
		return angle;
	}

	/**
	 * @brief Convert angles to degrees
	 */
	void ToDegree()
	{
		alpha = RadToDeg(alpha);
		beta = RadToDeg(beta);
		gamma = RadToDeg(gamma);
	}

	static double DegToRad(double deg)
	{
		return (deg*M_PI)/180;
	}

	static double RadToDeg(double rad)
	{
		return (rad*180)/M_PI;
	}
};

typedef Point3f Vector3f;
typedef Point3d Vector3d;

class Plane : public Polygon
{
public:
	Point3f normal;
};

class Geometry
{
public:
	static void DifferPolygons(const Polygon &subject, const Vector3f &subjNormal,
							   const Polygon &clip, const Vector3f &clipNormal,
							   const Vector3f &clipDir, PolygonArray &difference);

	static Point3f ProjectPointToPlane(const Point3f &point,
									   const Vector3f &direction,
									   const Vector3f &planeNormal);

	static bool IncidentBeamToFacet(Facet *facet, const Polygon &beamPol,
									bool isInside, const Vector3f &incDir,
									Polygon &intersection);

private:
	static bool ProjectPolygonToPlane(const Polygon &polygon, const Vector3f &dir,
									  const Point3f &normal, __m128 *_projection);

	static void RefineOutputPolygon(__m128 *_output_points, int outputSize,
									Polygon &polygon);
};
