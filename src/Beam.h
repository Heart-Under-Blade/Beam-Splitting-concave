#pragma once

#include "global.h"
#include "types.h"
#include "math/compl.hpp"
#include "JonesMatrix.h"
#include <vector> // DEB

class Beam
{
public:
	Beam();
	Beam(const Beam &other);

	void RotateSpherical(const Point3f &dir, const Point3f &polarBasis);
	void RotatePlane(const Point3f& newBasis); ///< rotate Jones matrix in case of beam splitting

	void AddVertex(const Point3f &vertex);
	void MulJMatrix(const Beam &other, const complex &coef1, const complex &coef2);
	void SetFormByOther(const Beam &other);

	Beam & operator = (const Beam &other);

public:
	Point3f direction;				///< direction of beam
	JonesMatrix JMatrix;			///< Jones matrix of beam
	Point3f e;						///< basis of polarization plane

	int facetId;					///< last reflected facet
	int level;						///< number of preview reflections
	bool isExternal;				///< beam state towards the particle (inside or outside)

	/// OPT: сделать др. вариант класса или трассировки
	/// и убрать эти параметры
	double opticalPath;				///< optical path of beam
	double D;						///< current position of phase front from Ax+By+Cz+D=0

	Point3f polygon[MAX_VERTEX_NUM];	///< array of beam vertices (shape)
	int size;							///< current vertex number

#ifdef _WRITE_TRACK
	std::vector<int> track;
#endif

private:
	void RotateJMatrix(const Point3f &newBasis);
	void GetSpherical(double &fi, double &teta) const;
	void Copy(const Beam &other);
};


