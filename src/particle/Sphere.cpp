#include "Sphere.h"

Sphere::Sphere(const complex &refrIndex, int nParallels, int nMeridians,
			   double radius)
{
	isConcave = false;
	Init((nParallels-1)*(nMeridians*2), refrIndex);

	Orientation angle;
	angle.beta = 0;
	angle.gamma = 0;

	double stepZen = M_PI/nParallels;
	double stepAz = 2*M_PI/nMeridians;

	int nFacets = 0;

	Point3f northPole(0, 0, radius);
	Point3f southPole(0, 0, -radius);

	std::vector<Point3f> prevParallelPoints;
	prevParallelPoints.resize(nMeridians);

	for (int iPar = 1; iPar <= nParallels; ++iPar)
	{
		angle.beta = iPar*stepZen;

		Point3f p1, p2, start;
		p1.coordinates[0] = radius*sin(angle.beta)*cos(angle.gamma);
		p1.coordinates[1] = radius*sin(angle.beta)*sin(angle.gamma);
		p1.coordinates[2] = radius*cos(angle.beta);

		for (int jMer = 1; jMer <= nMeridians; ++jMer)
		{
			angle.gamma = jMer*stepAz;

			if (iPar == 1) // pole
			{
				p2.coordinates[0] = radius*sin(angle.beta)*cos(angle.gamma);
				p2.coordinates[1] = radius*sin(angle.beta)*sin(angle.gamma);
				p2.coordinates[2] = radius*cos(angle.beta);

				defaultFacets[nFacets].AddVertex(northPole);
				defaultFacets[nFacets].AddVertex(p1);
				defaultFacets[nFacets].AddVertex(p2);
				++nFacets;

				prevParallelPoints[jMer-1] = p1;
			}
			else if (iPar == nParallels) // pole
			{
				defaultFacets[nFacets].AddVertex(southPole);
				defaultFacets[nFacets].AddVertex((jMer != nMeridians) ? prevParallelPoints[jMer] : start);
				defaultFacets[nFacets].AddVertex(prevParallelPoints[jMer-1]);
				++nFacets;
			}
			else
			{
				p2.coordinates[0] = radius*sin(angle.beta)*cos(angle.gamma);
				p2.coordinates[1] = radius*sin(angle.beta)*sin(angle.gamma);
				p2.coordinates[2] = radius*cos(angle.beta);

				defaultFacets[nFacets].AddVertex(prevParallelPoints[jMer-1]);
				defaultFacets[nFacets].AddVertex(p1);
				defaultFacets[nFacets].AddVertex(p2);
				++nFacets;

				defaultFacets[nFacets].AddVertex((jMer != nMeridians) ? prevParallelPoints[jMer] : start);
				defaultFacets[nFacets].AddVertex(prevParallelPoints[jMer-1]);
				defaultFacets[nFacets].AddVertex(p2);
				++nFacets;

				prevParallelPoints[jMer-1] = p1;
			}

			p1 = p2;
		}

		start = prevParallelPoints[0];
	}

	SetSymmetry(M_PI, M_PI*2);
	SetDefaultNormals();
	Reset();
	SetDefaultCenters();
}
