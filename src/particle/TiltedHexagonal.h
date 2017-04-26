//#pragma once

//#include "Particle.h"

///**
// * @brief The Hexagon class
// * The prism particle with 6 number of side facets.
// */
//class TiltedHexagonal : public Particle
//{
//public:
//	TiltedHexagonal();
//	TiltedHexagonal(double m_radius, double m_halfHeight, const complex &m_refractionIndex,
//					double tiltAngle);

////	void Rotate(double beta, double gamma, double alpha) override;

//private:
//	double m_tiltAngle;

//	void TiltBaseFacets();

//protected:
//	static const int BASE_FACET_NUM = 2;		///< number of bases of hexagon
//	static const int SIDE_VERTEX_NUMBER = 4;	///< number of vertex of the each side facet
//	static const int BASE_VERTEX_NUM = 6;		///< number of side facets

//	Polygon m_originBases[BASE_FACET_NUM];

//protected:
//	void SetDefaultNormals() override;
//	void SetFacetParams() override;

//	void SetBaseFacets();
//	void SetSideFacets(Point3f *baseTop, Point3f *baseBottom,
//					   int startIndex, int endIndex);
////	void RotateBaseFacets(Point3f *baseTop, Point3f *baseBottom);
////	void SetSideNormals(int beginId);
//};

