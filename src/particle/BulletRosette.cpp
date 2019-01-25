#include "BulletRosette.h"
#include "Bullet.h"
#include "common.h"

BulletRosette::BulletRosette()
{
}

BulletRosette::BulletRosette(const complex &refrIndex, const Size &size,
							 double peakHeight)
	: Particle(8, refrIndex, true) // REF: 8?????
{
	SetSymmetry(M_PI/2, M_PI);
	isAggregated = true;

	std::vector<Particle> bullets;
	double halfHeight = size.height/2;

	{
		Bullet b(refrIndex, size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle3d(0, Angle3d::DegToRad(90), 0));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle3d(0, Angle3d::DegToRad(180), 0));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle3d(0, Angle3d::DegToRad(270), 0));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle3d(0, Angle3d::DegToRad(90), 0));
		b.Fix();
		b.Rotate(Angle3d(Angle3d::DegToRad(90), 0, 0));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle3d(0, Angle3d::DegToRad(90), 0));
		b.Fix();
		b.Rotate(Angle3d(Angle3d::DegToRad(270), 0, 0));
		b.Fix();
		bullets.push_back(b);
	}

	Concate(bullets);

	SetDefaultNormals();
	SetDefaultCenters();
	Reset();

	for (int i = 0; i < nElems; ++i)
	{
		elems[i].origin.isOverlayedIn = true;
		elems[i].origin.isOverlayedOut = true;
		elems[i].actual.isOverlayedIn = true;
		elems[i].actual.isOverlayedOut = true;
	}
}

void BulletRosette::GetParticalFacetIdRange(Facet *facet, int &begin, int &end) const
{
	int patN = facet->index/13;
	begin = patN*13;
	end = begin+13;
}
