#pragma once

#include "Particle.h"
#include "global.h"
#include "math/compl.hpp"
#include "JonesMatrix.h"
#include "float.h"

#ifdef _TRACK_ALLOW
//#include <vector>
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

	// REF: перенести в PhisBeam
	complex DiffractionIncline(const Point3d& pt, double lam) const; ///< calculate diffraction at the point /b pt
	//--------------------------

	// REF: рассмотреть схему, где у пучка будет много полигонов

public:
	Point3f direction;				///< direction of beam
	JonesMatrix J;					///< Jones matrix of beam
	Point3f e;						///< basis of polarization plane

	Polygon polygon;				///< array of beam vertices (shape)

	int facetID;					///< last reflected facet
	int level;						///< number of preview reflections
	Location location;				///< beam state towards the particle (inside or outside)

	// REF: перенести в PhisBeam
//	Point3d T, F, N;
	double opticalPath;				///< optical path of beam
	double D;						///< current position of phase front from Ax+By+Cz+D=0

#ifdef _TRACK_ALLOW
	long long int id = 0;
//	std::vector<int> track;
#endif

//	friend Point3d Proj(const Point3d& _r, const Point3d &pnt)
//	{
//		Point3d _Tx,  // условная горизонталь СК экрана в СК тела
//				_Ty;  // третья ось (условная вертикаль СК экрана)
//		const double tmp = sqrt(SQR(_r.x)+SQR(_r.y));

//		(fabs(_r.z)>1-DBL_EPSILON) ? (_Tx=Point3d(0,-_r.z,0), _Ty=Point3d(1,0,0))
//								   : (_Tx=Point3d(_r.y/tmp,-_r.x/tmp,0), _Ty=CrossProductD(_r,_Tx));

//		return Proj(_Tx, _Ty, _r, pnt);
//	}

//	friend Point3d Proj(const Point3d& Tx, const Point3d& Ty,
//						const Point3d& r,  const Point3d& pnt)
//	{
//		const  Point3d p_pr = pnt-(r*DotProductD(r,pnt)); // расчёт коор-т в СК наблюдателя
//		return Point3d(DotProductD(p_pr,Tx), DotProductD(p_pr,Ty), 0); //*/
//	}

private:
	void RotateJMatrix(const Point3f &newBasis);
	void GetSpherical(double &fi, double &teta) const;
	void Copy(const Beam &other);
};


