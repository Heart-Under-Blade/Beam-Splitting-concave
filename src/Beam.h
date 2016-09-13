#pragma once

#include "global.h"
#include "math/compl.hpp"
#include "JonesMatrix.h"

class Beam
{
public:
	Beam();
	Beam(const Beam &other);

	void RotateSpherical(const Point3f &dir, const Point3f &polarBasis);
	void RotatePlane(const Point3f& newBasis); ///< rotate Jones matrix in case of beam splitting

	double Square() const;
	Point3f Center() const;
	double CrossSection() const;
	void AddVertex(const Point3f &vertex);
	void MulJMatrix(const Beam &other, const complex &coef1, const complex &coef2);

	Beam & operator = (const Beam &other);

public:
	Point3f direction;				///< direction of the beam in 3D space
	Point3f e;						///< basis of polarization plane
	JonesMatrix JMatrix;				///< Jones matrix of the beam
	double D;						///< current position of phase front from Ax+By+Cz+D=0
	double opticalPath;				///< optical path of the beam

	Point3f shape[MAX_VERTEX_NUM];	///< beam's vertices
	int shapeSize;					///< current vertex number of shape

private:
	void RotateJMatrix(const Point3f &newBasis);
	void GetSpherical(double &fi, double &teta) const;
	void copy(const Beam &other);
};
