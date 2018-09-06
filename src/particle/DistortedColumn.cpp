#include "DistortedColumn.h"
#include "global.h"
#include <iostream>

DistortedColumn::DistortedColumn(const complex &refrIndex, const Size &size,
								 double angle)
	: Column(8, refrIndex, size, false)
{
	SetSymmetry(M_PI/2, 2*M_PI);
	SetFacetParams();

	SetBases(defaultFacets[0], defaultFacets[7]);
	DistortBases(angle);
	SetSides(defaultFacets[0], defaultFacets[7]);

	SetDefaultNormals();
	SetDefaultCenters();
	Reset();
}

void DistortedColumn::DistortBases(double angle)
{
    double tilt = Angle::DegToRad(15);
    double tanA = tan(Angle::DegToRad(angle));
	double k = m_size.diameter/2 * tanA;

	double h[6];
	h[0] = k * cos(M_PI/3.0-tilt);
	h[1] = k * cos(2.0*M_PI/3.0-tilt);
	h[2] = k * cos(M_PI-tilt);
	h[3] = k * cos(4.0*M_PI/3.0-tilt);
	h[4] = k * cos(5.0*M_PI/3.0-tilt);
	h[5] = k * cos(tilt);

	int endPointIndex = BASE_VERTEX_NUM-1;

    for (size_t i = 0; i < nFacets; ++i)
	{
		defaultFacets[0].arr[i].cz += h[i];
		defaultFacets[7].arr[endPointIndex-i].cz += h[i];
	}
}
