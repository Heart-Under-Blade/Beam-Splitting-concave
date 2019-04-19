#pragma once

#include "Column.h"

class Bullet : public Column
{
public:
	Bullet();
	Bullet(const Size &size, double peakHeight);

private:
	void SetBaseFacet(Facet &facet);
	void SetPeakFacets(int start, int end, const Point3f *baseFacet,
					   const Point3f &peakPoint);

protected:
	void SetFacetParams() override;
};
