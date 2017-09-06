#include "BulletRosette.h"
#include "Bullet.h"
#include "global.h"

BulletRosette::BulletRosette()
{

}

BulletRosette::BulletRosette(const complex &refrIndex, double diameter,
							 double height, double peakHeight)
{
	double size = height*2 + diameter;
	Init(8, refrIndex, size);

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
		b.Rotate(DegToRad(90), 0, 0);
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(DegToRad(180), 0, 0);
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(DegToRad(270), 0, 0);
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -(halfHeight + peakHeight + 1));
		b.Fix();
		b.Rotate(DegToRad(90), 0, 0);
		b.Fix();
		b.Rotate(0, 0, DegToRad(90));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
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
	Reset();
	SetDefaultCenters();
}
