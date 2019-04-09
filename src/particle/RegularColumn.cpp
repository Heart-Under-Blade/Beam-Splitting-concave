#include "RegularColumn.h"

RegularColumn::RegularColumn(const Size &size)
	: Column(8, size, false)
{
	SetSymmetry(M_PI/2, M_PI/3);
	SetFacetParams();

	SetBases(elems[0].original, elems[7].original);
	SetSides(elems[0].original, elems[7].original);

	SetDefaultNormals();
	ResetPosition();
	SetDefaultCenters();
}
