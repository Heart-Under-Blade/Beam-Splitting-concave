#pragma once

#include "Column.h"

class Bullet : public Column
{
public:
	Bullet();
	Bullet(const complex &refrIndex, double diameter, double height, double peakHeight);

private:
	void SetBaseFacet(Facet &facet);
	void SetPeakFacets(int start, int end, const Point3f *baseFacet,
					   const Point3f &peakPoint);

protected:
	void SetFacetParams() override;
};
