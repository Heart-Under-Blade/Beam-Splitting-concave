#pragma once

#include "global.h"
#include "math/compl.hpp"
#include "JonesMatrix.h"
#include "float.h"
#include "BigInteger.hh"
#include "geometry_lib.h"
#include "Polygon.h"

class Light
{
public:
	Point3f direction;
	Point3f polarizationBasis;
};

class Track
{
public:
	BigInteger id = 0;
	int locations;		///< each bit of variable represents location of beam after an r/r act from left to right
						///< "0" when beam location is "inside" and "1" if it's "outside"

	Location GetLocationByActNumber(int nActs) const
	{
		int mask = 1;
		mask <<= nActs;
		return (locations & mask) ? Location::Out : Location::In;
	}
};

class Beam : public Polygon, public Light, public Track
{
public:
	Beam();
	Beam(const Beam &other);
	Beam(const Polygon &other);
	Beam(Beam &&other);

	void RotateSpherical(const Vector3f &dir, const Vector3f &polarBasis);

	void SetPolygon(const Polygon &other);
	void SetLight(const Vector3f &dir, const Vector3f &polarBasis);
	void SetLight(const Light &other);
	void AddOpticalPath(double path);
	void CopyTrack(const Track &other);

	Beam & operator = (const Beam &other);
	Beam & operator = (const Polygon &other);
	Beam & operator = (const Light &other);
	Beam & operator = (Beam &&other);

	void SetTracingParams(int facetId, int actN, Location location);

	void MultiplyJonesMatrix(const complex &c1, const complex &c2);
	void RotateJMatrix(const Vector3f &newBasis);

	// REF: перенести в PhisBeam
	complex DiffractionIncline(const Point3d& pt, double wavelength) const; ///< calculate diffraction at the point /b pt
	//--------------------------

	friend std::ostream & operator << (std::ostream &os, const Beam &beam);

	// REF: рассмотреть схему, где у пучка будет много полигонов

public:
	Matrix2x2c J;		///< Jones matrix of beam

	int nActs;			///< number of preview reflections
	int lastFacetId;	///< last reflected facet id
	Location location; // REF: заменить на 'bool isInside'			///< beam state towards the particle (inside or outside)

	// REF: перенести в PhisBeam
	double opticalPath;	///< optical path of beam
	double front;		///< current position of phase front from Ax+By+Cz+D=0 (where D is front)

#ifdef _DEBUG // DEB
	std::vector<Point3f> dirs;
	std::vector<double> ops;
#endif

#ifdef _TRACK_ALLOW
	BigInteger trackId = 0;
//	long long trackId = 0;
#endif

private:
	void GetSpherical(double &fi, double &teta) const;
	void Copy(const Beam &other);
	void SetDefault(Beam &other);
};
