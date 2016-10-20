#pragma once

#include "Hexagonal.h"

/**
 * @brief The ConcaveHexagonal class
 * The prism particle with 6 number of side facets and 2 cavities on the base facets.
 */
class ConcaveHexagonal : public Hexagonal
{
public:
	ConcaveHexagonal(double m_radius, double halfHeight, const complex &refractionIndex,
					 double cavityDept);

	void Rotate(double beta, double gamma, double alpha) override;
/**
 * OPT: сделать список граней, пучок, вышедший из которых,
 * не встретит затеняющих граней.
 * + добавить такую проверку в трассировку, чтобы избежать лишних операций
  */
private:
	double m_cavityDept;

	struct CavityPoints
	{
		Point3f top;
		Point3f bottom;
	}
	m_originCavities;

	const int CAVITY_FACET_VERTEX_NUM = 3;	///< number of vertex of the each facet into the cavity

protected:
	void SetFacetParams() override;
	void SetOriginNormals() override;

private:
	void SetCavities(Point3f *baseTop, Point3f *baseBottom, const CavityPoints &cavities);

	void SetCavityFacets(int start, int end, Point3f *baseFacet,
						 const Point3f &cavityPoint);
	void SetOriginCavityPoints();
	void SetBaseNormals();
	void SetSideNormals();
};
