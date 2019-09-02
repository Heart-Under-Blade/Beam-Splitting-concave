#pragma once

#include "Facet.h"
#include "compl.hpp"
#include "JonesMatrix.h"
#include "float.h"
#include "BigInteger.hh"
#include "geometry_lib.h"
#include "Polygon.h"
#include "TrackTree.h"

#ifdef _DEBUG // DEB
typedef long long IdType;
#else
typedef BigInteger IdType;
#endif

class Light
{
public:
	Point3f direction;
	Point3f polarizationBasis;
};

class Track
{
public:
	IdType id; ///< Unique id of beam calculated by facet ids

	// OPT бесполезно для выпуклых частиц (там всегда пучок внутри)
	int locations;	///< Every bit of the variable represents location of beam
					///< after an r/r act from the right to the left
					///< "0" when beam location is "inside" and "1" if it's "outside"

	int actNo; ///< Current r/r act number
	Facet *facet; ///< Last incident facet of the Particle
	TrackNode *node;

	Track()
	{
		id = 0;
		actNo = -1;
		locations = 0;
	}

	void Update(Facet *face)
	{
		facet = face;
		++actNo;
	}

	/**
	 * @brief Shows where the Beam is located towards the Particle inside or outside
	 * @param iAct r/r act number
	 * @return true if beam location is "inside" on r/r act, otherwise returns false
	 */
	bool IsInsideOnAct(int iAct) const
	{
		int mask = 1;
		mask <<= iAct;
		return !(locations & mask);
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

	void SetMatrix(const Matrix2x2c &matrix);
	void SetPolygon(const Polygon &other);
	void SetLight(const Vector3f &dir, const Vector3f &polarBasis);
	void SetLight(const Light &other);
	void AddOpticalPath(double path);
	void CopyTrack(const Track &other);

	Beam & operator = (const Beam &other);
	Beam & operator = (const Polygon &other);
	Beam & operator = (const Light &other);
	Beam & operator = (Beam &&other);

	void MultiplyJonesMatrix(const complex &c1, const complex &c2);
	void RotateJMatrix(const Vector3f &newBasis);

	//--------------------------

	friend std::ostream & operator << (std::ostream &os, const Beam &beam);

	// REF: рассмотреть схему, где у пучка будет много полигонов

public:
	Matrix2x2c J;		///< Jones matrix of beam

	int nActs;			///< number of preview reflections
	int lastFacetId;	///< last reflected facet id
	bool isInside; // REF: заменить на 'bool isInside'			///< beam state towards the particle (inside or outside)

	// REF: перенести в PhisBeam
	double opticalPath;	///< optical path of beam
	double front;		///< current position of phase front from Ax+By+Cz+D=0 (where D is front)

#ifndef _DEBUG // DEB
//	std::vector<Polygon> pols;
//	std::vector<Point3f> dirs;
	std::vector<double> ops;
#endif

private:
	void GetSpherical(double &fi, double &teta) const;
	void Copy(const Beam &other);
	void SetDefault(Beam &other);
};
