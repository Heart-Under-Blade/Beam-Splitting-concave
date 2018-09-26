#pragma once

#include "global.h"
#include "compl.hpp"
#include "JonesMatrix.h"
#include "float.h"
#include "BigInteger.hh"
#include "geometry_lib.h"
#include "Facet.h"
#include "Tracks.h"

class Light
{
public:
	Point3f direction;
	Point3f polarizationBasis;
};

class Track
{
public:
	IdType id = 0;
	int locations;		///< each bit of variable represents location of beam after an r/r act from left to right
						///< "0" when beam location is "inside" and "1" if it's "outside"

	bool IsInsideOnAct(int act) const
	{
		int mask = 1;
		mask <<= act;
		return !(locations & mask);
	}

	Track & operator = (const Track &other)
	{
		if (this != &other) // OPT: попробовать убрать это уловие для ускорения
		{
			id = other.id;
			locations = other.locations;
		}

		return *this;
	}
};

class Beam : public Polygon, public Light, public Track
{
public:
	Beam();
	Beam(const Beam &other);
	Beam(const Polygon &other);
	Beam(Beam &&other);

	Vector3f RotateSpherical(const Vector3f &dir, const Vector3f &polarBasis);

	void SetPolygon(const Polygon &other);
	void SetLight(const Vector3f &dir, const Vector3f &polarBasis);
	void SetLight(const Light &other);
	void Clear();
	void AddOpticalPath(double path);
	void CopyTrack(const Track &other);

	Beam & operator = (const Beam &other);
	Beam & operator = (const Polygon &other);
	Beam & operator = (const Light &other);
	Beam & operator = (Beam &&other);

	void SetTracingParams(Facet *fac, int actN, bool isIn);

	void MultiplyJonesMatrix(const complex &f1, const complex &f2);
	void RotateJMatrix(const Vector3f &newBasis);

	// REF: перенести в PhisBeam
	complex DiffractionIncline(const Point3d& pt, double wavelength) const; ///< calculate diffraction at the point /b pt
	//--------------------------

	friend std::ostream & operator << (std::ostream &os, const Beam &beam);

	// REF: рассмотреть схему, где у пучка будет много полигонов

public:
	Matrix2x2c J;		///< Jones matrix of beam
	int act;			///< number of preview reflections
	Facet *facet;		///< last incident facet

	bool isInside; 		///< beam state towards the particle (inside or outside)

	// REF: перенести в PhisBeam
	double opticalPath;	///< optical path of beam
	double front;		///< current position of phase front from Ax+By+Cz+D=0 (where D is front)

#ifdef _DEBUG // DEB
	std::vector<Polygon> pols;
	std::vector<Point3f> dirs;
	std::vector<double> ops;
#endif

private:
	void GetSpherical(double &fi, double &teta) const;
	void Copy(const Beam &other);
	void SetDefault(Beam &other);
};
