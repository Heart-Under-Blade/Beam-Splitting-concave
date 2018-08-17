#include "DistortedHexagonal.h"
#include "global.h"
#include <iostream>

DistortedHexagonal::DistortedHexagonal(const complex &refrIndex, double diameter, double height,
									   double angle)
{
	isConcave = false;
	SetSize(diameter, height);
	Init(8, refrIndex); // REF: перенести в конструктор Hexagonal и вызывать его тут

	SetSymmetry(M_PI/2, 2*M_PI);
	SetFacetParams();

	SetBases(defaultFacets[0], defaultFacets[7]);
	DistortBases(angle);
	SetSides(defaultFacets[0], defaultFacets[7]);

	SetDefaultNormals();
	SetDefaultCenters();
	Reset();
}

void DistortedHexagonal::DistortBases(double angle)
{
	double tilt = DegToRad(15);
	double tanA = tan(DegToRad(angle));
	double k = m_diameter/2 * tanA;

	double h[6];
	h[0] = k * cos(5.0*M_PI/3.0-tilt);
	h[1] = k * cos(tilt);
	h[2] = k * cos(M_PI/3.0-tilt);
	h[3] = k * cos(2.0*M_PI/3.0-tilt);
	h[4] = k * cos(M_PI-tilt);
	h[5] = k * cos(4.0*M_PI/3.0-tilt);

	int endPointIndex = BASE_VERTEX_NUM-1;

	for (int i = 0; i < nFacets; ++i)
	{
		defaultFacets[0].arr[i] += h[i];
		defaultFacets[7].arr[endPointIndex-i] += h[i];
	}
}
