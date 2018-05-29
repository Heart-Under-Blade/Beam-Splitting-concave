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

class Beam : public Polygon, public Light
{
public:
	Beam();
	Beam(const Beam &other);
	Beam(const Polygon &other);
	Beam(Beam &&other);

	void RotateSpherical(const Vector3f &dir, const Vector3f &polarBasis);
	void RotatePlane(const Point3f& newBasis); ///< rotate Jones matrix in case of beam splitting

	Location GetLocationByActNumber(int act) const;

	void AddVertex(const Point3f &vertex);
	void SetPolygon(const Polygon &other);
	void SetLight(const Point3f &dir, const Point3f &polarBasis);
	void ComputeFrontPosition();

	Beam & operator = (const Beam &other);
	Beam & operator = (const Polygon &other);
	Beam & operator = (const Light &other);
	Beam & operator = (Beam &&other);

	void SetTracingParams(int facetID, int act, Location location);

	void SetJonesMatrix(const Beam &other, const complex &c1, const complex &c2);

	// REF: перенести в PhisBeam
	complex DiffractionIncline(const Point3d& pt, double wavelength) const; ///< calculate diffraction at the point /b pt
	//--------------------------

	friend std::ostream & operator << (std::ostream &os, const Beam &beam);

	// REF: рассмотреть схему, где у пучка будет много полигонов

public:
	Matrix2x2c J;					///< Jones matrix of beam

	int lastFacetId;				///< last reflected facet
	int act;						///< number of preview reflections
	Location location;				///< beam state towards the particle (inside or outside)

	// REF: перенести в PhisBeam
	double opticalPath;				///< optical path of beam
	double frontPosition;			///< current position of phase front from Ax+By+Cz+D=0

#ifdef _TRACK_ALLOW
	BigInteger trackId = 0;
//	std::vector<int> track;
#endif

#ifdef _DEBUG // DEB
	std::vector<Point3f> dirs;
	std::vector<double> ops;
#endif
	int locations;					///< each bit of variable represents location of beam after an r/r act from left to right
									///< "0" when beam location is "inside" and "1" if it's "outside"

private:
	void RotateJMatrix(const Point3f &newBasis);
	void GetSpherical(double &fi, double &teta) const;
	void Copy(const Beam &other);
	void SetDefault(Beam &other);
};
