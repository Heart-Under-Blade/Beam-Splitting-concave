#include "BulletRosette.h"
#include "Bullet.h"
#include "global.h"

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
		b.Rotate(Angle(0, Angle::DegToRad(90), 0));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle(0, Angle::DegToRad(180), 0));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle(0, Angle::DegToRad(270), 0));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, size, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(Angle(0, Angle::DegToRad(90), 0));
		b.Fix();
		b.Rotate(Angle(Angle::DegToRad(90), 0, 0));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, size, peakHeight);
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

	for (size_t i = 0; i < nElems; ++i)
	{
		elems[i].origin.isOverlayedIn = false;
		elems[i].origin.isOverlayedOut = false;
		elems[i].actual.isOverlayedIn = false;
		elems[i].actual.isOverlayedOut = false;
	}
}

void BulletRosette::GetParticalFacetIdRangeByFacetId(int id, int &begin, int &end) const
{
	int patN = id/13;
	begin = patN*13;
	end = begin+13;
}
