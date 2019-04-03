#pragma once

#include "Column.h"

/**
 * @brief The ConcaveHexagonal class
 * The prism particle with 6 number of side facets and 2 cavities on the base facets.
 */
class HollowColumn : public Column
{
public:
	HollowColumn(const complex &refrIndex, const Size &size,
				 double cavityAngle);

private:
	double m_cavityDept;

	struct CavityPoints
	{
		Point3f top;
		Point3f bottom;
	}
	m_defaultStateCavities;

	const int CAVITY_FACET_VERTEX_NUM = 3;	///< number of vertex of the each facet into the cavity

protected:
	void SetFacetParams() override;

private:
	void SetCavities(Facet &baseTop, Facet &baseBottom, const CavityPoints &cavities);

	void SetCavityFacets(int start, int end, const Point3f *baseFacet,
						 const Point3f &cavityPoint);
	void SetOriginCavityPoints();

};
