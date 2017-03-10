#pragma once

#include "Hexagonal.h"

/**
 * @brief The ConcaveHexagonal class
 * The prism particle with 6 number of side facets and 2 cavities on the base facets.
 */
class ConcaveHexagonal : public Hexagonal
{
public:
	ConcaveHexagonal(double radius, double halfHeight, const complex &refractionIndex,
					 double cavityDept);

	void Rotate(double beta, double gamma, double alpha) override;

private:
	double m_cavityDept;

	struct CavityPoints
	{
		Point3f top;
		Point3f bottom;
	}
	m_defaultCavities;

	const int CAVITY_FACET_VERTEX_NUM = 3;	///< number of vertex of the each facet into the cavity

protected:
	void SetFacetParams() override;
	void SetDefaultNormals() override;

private:
	void SetCenters();
	void SetCavities(const CavityPoints &cavities);

	void SetCavityFacets(int start, int end, const Point3f *baseFacet,
						 const Point3f &cavityPoint);
	void SetOriginCavityPoints();
	void SetBaseNormals();
	void RotateCavities();
};
