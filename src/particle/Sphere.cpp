#include "Sphere.h"

Sphere::Sphere(int nParallels, int nMeridians, double radius)
	: Particle((nParallels-1)*(nMeridians*2), false)
{

	Orientation angle;
	angle.zenith = 0;
	angle.azimuth = 0;

	double stepZen = M_PI/nParallels;
	double stepAz = 2*M_PI/nMeridians;

	int nFacets = 0;

	Point3f northPole(0, 0, radius);
	Point3f southPole(0, 0, -radius);

	std::vector<Point3f> prevParallelPoints;
	prevParallelPoints.resize(nMeridians);

	for (int iPar = 1; iPar <= nParallels; ++iPar)
	{
		angle.zenith = iPar*stepZen;

		Point3f p1, p2, start;
		p1.coordinates[0] = radius*sin(angle.zenith)*cos(angle.azimuth);
		p1.coordinates[1] = radius*sin(angle.zenith)*sin(angle.azimuth);
		p1.coordinates[2] = radius*cos(angle.zenith);

		for (int jMer = 1; jMer <= nMeridians; ++jMer)
		{
			angle.azimuth = jMer*stepAz;

			if (iPar == 1) // pole
			{
				p2.coordinates[0] = radius*sin(angle.zenith)*cos(angle.azimuth);
				p2.coordinates[1] = radius*sin(angle.zenith)*sin(angle.azimuth);
				p2.coordinates[2] = radius*cos(angle.zenith);

				elems[nFacets].original.AddVertex(northPole);
				elems[nFacets].original.AddVertex(p1);
				elems[nFacets].original.AddVertex(p2);
				++nFacets;

				prevParallelPoints[jMer-1] = p1;
			}
			else if (iPar == nParallels) // pole
			{
				elems[nFacets].original.AddVertex(southPole);
				elems[nFacets].original.AddVertex((jMer != nMeridians) ? prevParallelPoints[jMer] : start);
				elems[nFacets].original.AddVertex(prevParallelPoints[jMer-1]);
				++nFacets;
			}
			else
			{
				p2.coordinates[0] = radius*sin(angle.zenith)*cos(angle.azimuth);
				p2.coordinates[1] = radius*sin(angle.zenith)*sin(angle.azimuth);
				p2.coordinates[2] = radius*cos(angle.zenith);

				elems[nFacets].original.AddVertex(prevParallelPoints[jMer-1]);
				elems[nFacets].original.AddVertex(p1);
				elems[nFacets].original.AddVertex(p2);
				++nFacets;

				elems[nFacets].original.AddVertex((jMer != nMeridians) ? prevParallelPoints[jMer] : start);
				elems[nFacets].original.AddVertex(prevParallelPoints[jMer-1]);
				elems[nFacets].original.AddVertex(p2);
				++nFacets;

				prevParallelPoints[jMer-1] = p1;
			}

			p1 = p2;
		}

		start = prevParallelPoints[0];
	}

	SetSymmetry(M_PI, M_PI*2);
	SetDefaultNormals();
	ResetPosition();
	SetDefaultCenters();
}
