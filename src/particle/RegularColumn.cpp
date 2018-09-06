#include "RegularColumn.h"

RegularColumn::RegularColumn(const complex &refrIndex, const Size &size)
	: Column(8, refrIndex, size, false)
{
	SetSymmetry(M_PI/2, M_PI/3);
	SetFacetParams();

	SetBases(defaultFacets[0], defaultFacets[7]);
	SetSides(defaultFacets[0], defaultFacets[7]);

	SetDefaultNormals();
	Reset();
	SetDefaultCenters();
}
