#include "BulletRosette.h"
#include "Bullet.h"
#include "global.h"

BulletRosette::BulletRosette()
{
	SetSymmetry(M_PI/2, M_PI);
}

<<<<<<< HEAD
BulletRosette::BulletRosette(const complex &refrIndex, double diameter,
							 double height, double peakHeight)
	: BulletRosette()
=======
BulletRosette::BulletRosette(const Size &size, double peakHeight)
	: Particle(8, true) // REF: 8?????
>>>>>>> origin/refactor
{
	isConcave = true;
	isAggregated = true;

	Init(8, refrIndex);

	std::vector<Particle> bullets;
	double halfHeight = height/2;

	{
<<<<<<< HEAD
		Bullet b(refrIndex, diameter, height, peakHeight);
=======
		Bullet b(size, peakHeight);
>>>>>>> origin/refactor
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		bullets.push_back(b);
	}

	{
<<<<<<< HEAD
		Bullet b(refrIndex, diameter, height, peakHeight);
=======
		Bullet b(size, peakHeight);
>>>>>>> origin/refactor
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(DegToRad(90), 0, 0);
		b.Fix();
		bullets.push_back(b);
	}

	{
<<<<<<< HEAD
		Bullet b(refrIndex, diameter, height, peakHeight);
=======
		Bullet b(size, peakHeight);
>>>>>>> origin/refactor
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(DegToRad(180), 0, 0);
		b.Fix();
		bullets.push_back(b);
	}

	{
<<<<<<< HEAD
		Bullet b(refrIndex, diameter, height, peakHeight);
=======
		Bullet b(size, peakHeight);
>>>>>>> origin/refactor
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(DegToRad(270), 0, 0);
		b.Fix();
		bullets.push_back(b);
	}

	{
<<<<<<< HEAD
		Bullet b(refrIndex, diameter, height, peakHeight);
=======
		Bullet b(size, peakHeight);
>>>>>>> origin/refactor
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(DegToRad(90), 0, 0);
		b.Fix();
		b.Rotate(0, 0, DegToRad(90));
		b.Fix();
		bullets.push_back(b);
	}

	{
<<<<<<< HEAD
		Bullet b(refrIndex, diameter, height, peakHeight);
=======
		Bullet b(size, peakHeight);
>>>>>>> origin/refactor
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(DegToRad(90), 0, 0);
		b.Fix();
		b.Rotate(0, 0, DegToRad(270));
		b.Fix();
		bullets.push_back(b);
	}

	Concate(bullets);

	SetDefaultNormals();
	SetDefaultCenters();
	Reset();

	for (int i = 0; i < nFacets; ++i)
	{
		defaultFacets[i].isVisibleIn = false;
		defaultFacets[i].isVisibleOut = false;
		facets[i].isVisibleIn = false;
		facets[i].isVisibleOut = false;
	}
}

<<<<<<< HEAD
void BulletRosette::GetParticalFacetIdRangeByFacetId(int id, int &begin, int &end) const
{
	int patN = id/13;
	begin = patN*13;
	end = begin+13;
=======
void BulletRosette::GetPartByFacet(Facet *facet, Array<Facet*> &facets)
{
	int patN = facet->index/13;
	int begin = patN*13;
	int end = begin+13;

	GetFacets(begin, end, facets);
>>>>>>> origin/refactor
}
