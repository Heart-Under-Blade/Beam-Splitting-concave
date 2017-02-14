#pragma once

#include "Particle.h"
#include "global.h"
#include "math/compl.hpp"
#include "JonesMatrix.h"

#ifdef _TRACK_ALLOW
#include <vector>
#endif

class Beam
{
public:
	Beam();
	Beam(const Beam &other);
	Beam(Beam &&other);

	void RotateSpherical(const Point3f &dir, const Point3f &polarBasis);
	void RotatePlane(const Point3f& newBasis); ///< rotate Jones matrix in case of beam splitting

	void AddVertex(const Point3f &vertex);
	void SetPolygon(const Polygon &other);

	Beam & operator = (const Beam &other);
	Beam & operator = (Beam &&other);

	// REF: рассмотреть схему, где у пучка будет много полигонов

public:
	Point3f direction;				///< direction of beam
	JonesMatrix JMatrix;			///< Jones matrix of beam
	Point3f e;						///< basis of polarization plane

	Polygon polygon;				///< array of beam vertices (shape)

	int facetID;					///< last reflected facet
	int level;						///< number of preview reflections
	Location location;				///< beam state towards the particle (inside or outside)

	/// OPT: сделать др. вариант класса или трассировки
	/// и убрать эти параметры
	double opticalPath;				///< optical path of beam
	double D;						///< current position of phase front from Ax+By+Cz+D=0

#ifdef _TRACK_ALLOW
	std::vector<int> track;
#endif

private:
	void RotateJMatrix(const Point3f &newBasis);
	void GetSpherical(double &fi, double &teta) const;
	void Copy(const Beam &other);
};


