#include "BulletRosette.h"
#include "Bullet.h"
#include "common.h"
#include "geometry_lib.h"

BulletRosette::BulletRosette()
{
}

BulletRosette::BulletRosette(const Size &size, double peakHeight)
	: Particle(8, true) // REF: 8?????
{
	SetSymmetry(M_PI/2, M_PI);
	isAggregated = true;

	std::vector<Particle> bullets;
	double halfHeight = size.height/2;

	{
		Bullet b(size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.CommitState();
		bullets.push_back(b);
	}

	{
		Bullet b(size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.CommitState();
		b.Rotate(Orientation(0, Orientation::DegToRad(90)));
		b.CommitState();
		bullets.push_back(b);
	}

	{
		Bullet b(size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.CommitState();
		b.Rotate(Orientation(0, Orientation::DegToRad(180)));
		b.CommitState();
		bullets.push_back(b);
	}

	{
		Bullet b(size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.CommitState();
		b.Rotate(Orientation(0, Orientation::DegToRad(270)));
		b.CommitState();
		bullets.push_back(b);
	}

	{
		Bullet b(size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.CommitState();
		b.Rotate(Orientation(0, Orientation::DegToRad(90)));
		b.CommitState();
		b.Rotate(Orientation(Orientation::DegToRad(90), 0));
		b.CommitState();
		bullets.push_back(b);
	}

	{
		Bullet b(size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.CommitState();
		b.Rotate(Orientation(0, Orientation::DegToRad(90)));
		b.CommitState();
		b.Rotate(Orientation(Orientation::DegToRad(270), 0));
		b.CommitState();
		bullets.push_back(b);
	}

	Concatenate(bullets);

	SetDefaultNormals();
	SetDefaultCenters();
	ResetPosition();

	for (int i = 0; i < nElems; ++i)
	{
		elems[i].original.isOverlayedIn = true;
		elems[i].original.isOverlayedOut = true;
		elems[i].actual.isOverlayedIn = true;
		elems[i].actual.isOverlayedOut = true;
	}
}

void BulletRosette::GetPartByFacet(Facet *facet, Array<Facet*> &facets)
{
	int patN = facet->index/13;
	int begin = patN*13;
	int end = begin+13;

	GetFacets(begin, end, facets);
}
