#include "BulletRosette.h"
#include "Bullet.h"
#include "global.h"

BulletRosette::BulletRosette()
{
	SetSymmetry(M_PI/2, M_PI);
}

BulletRosette::BulletRosette(const complex &refrIndex, double diameter,
							 double height, double peakHeight)
	: BulletRosette()
{
	isConcave = true;
	isAggregated = true;

	Init(8, refrIndex);

	std::vector<Particle> bullets;
	double halfHeight = height/2;

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle(0, Angle::DegToRad(90), 0));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle(0, Angle::DegToRad(180), 0));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle(0, Angle::DegToRad(270), 0));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle(0, Angle::DegToRad(90), 0));
		b.Fix();
		b.Rotate(Angle(Angle::DegToRad(90), 0, 0));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle(0, Angle::DegToRad(90), 0));
		b.Fix();
		b.Rotate(Angle(Angle::DegToRad(270), 0, 0));
		b.Fix();
		bullets.push_back(b);
	}

	Concate(bullets);

	SetDefaultNormals();
	SetDefaultCenters();
	Reset();

	for (size_t i = 0; i < nFacets; ++i)
	{
		defaultFacets[i].isVisibleIn = false;
		defaultFacets[i].isVisibleOut = false;
		facets[i].isVisibleIn = false;
		facets[i].isVisibleOut = false;
	}
}

void BulletRosette::GetParticalFacetIdRangeByFacetId(int id, int &begin, int &end) const
{
	int patN = id/13;
	begin = patN*13;
	end = begin+13;
}
