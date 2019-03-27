#include "RegularColumn.h"

RegularColumn::RegularColumn(const complex &refrIndex, const Size &size)
	: Column(8, refrIndex, size, false)
{
	SetSymmetry(M_PI/2, M_PI/3);
	SetFacetParams();

	SetBases(elems[0].origin, elems[7].origin);
	SetSides(elems[0].origin, elems[7].origin);

	SetDefaultNormals();
	Reset();
	SetDefaultCenters();
}
