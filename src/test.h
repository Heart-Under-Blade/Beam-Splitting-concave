/**
  File for various tests
*/

#pragma once

#include "global.h"
#include <fstream>
#include <iostream>

#include "Hexagonal.h"
#include "HexagonalAggregate.h"
#include "TiltedHexagonal.h"
#include "ConcaveHexagonal.h"
#include "Intersection.h"
#include "Scattering.h"

void outputParticle(const Particle &particle)
{
	for (int i = 0; i < particle.nFacets; ++i)
	{
		std::cout << i << ": ";
		for (int j = 0; j < particle.facets[i].size; ++j)
		{
			Point3f p = particle.facets[i].arr[j];
			std::cout << "("
					  << p.point[0] << ", "
					  << p.point[1] << ", "
					  << p.point[2]
					  << "), ";
		}

		std::cout << std::endl;
	}

	std::cout << std::endl << "Normals" << std::endl << std::endl;

	for (int i = 0; i < particle.nFacets; ++i)
	{
		std::cout << i << ": ";
		std::cout << "("
				  << particle.facets[i].in_normal.cx << ", "
				  << particle.facets[i].in_normal.cy << ", "
				  << particle.facets[i].in_normal.cz
				  << "), ";
		std::cout << std::endl;
	}
}

void toFile(const Particle &particle)
{
	std::ofstream M("particle.dat", std::ios::out);

	for (int i = 0; i < particle.nFacets; ++i)
	{
		for (int j = 0; j < particle.facets[i].size; ++j)
		{
			Point3f p = particle.facets[i].arr[j];
			M << p.point[0] << ' '
							<< p.point[1] << ' '
							<< p.point[2];
			M << std::endl;
		}

		M << std::endl << std::endl;;
	}

	M.close();
}

void testHexagonBuilding()
{
//	Particle hex = Hexagonal(40, 100, 1.31);
//	toFile(hex);
//	outputParticle(hex);
}

void testConcaveHexagonRot()
{
//	Particle *hex = new ConcaveHexagonal(20, 100, 1.31, 10);
//	double beta = ((19 + 0.5)*M_PI)/(2.0*100);
//	double gamma = ((36 + 0.5)*M_PI)/(3.0*101);
//	hex->Rotate(beta, gamma, 0);
//	toFile(*hex);
//	outputParticle(*hex);
}

void testHexagonRotate()
{
//	Particle *hex = new Hexagonal(40, 100, 1.31);
//	double beta = 19*M_PI/180;
//	double gamma = 90*M_PI/180;
//	double beta = ((19 + 0.5)*M_PI)/(2.0*100);
//	double gamma = ((36 + 0.5)*M_PI)/(3.0*101);
//	hex->Rotate(beta, gamma, 0);

//	toFile(*hex);
//	outputParticle(*hex);
}


//void testTiltHexagonBuild()
//{
//	Particle *hex = new TiltedHexagonal(34.8, 50, 1.31, 0.25);
//	double beta = ((19 + 0.5)*M_PI)/(2.0*100);
//	double gamma = ((36 + 0.5)*M_PI)/(3.0*101);
////	double beta = 90*M_PI/180;
////	double gamma = 0*M_PI/180;
////	double beta = ((19 + 0.5)*M_PI)/(2.0*100);
////	double gamma = ((36 + 0.5)*M_PI)/(3.0*101);
////	hex->Rotate(beta, gamma, 0);

//	toFile(*hex);
//	outputParticle(*hex);
//}

void testHexagonalAggregateBuild()
{
	Particle *hex = new HexagonalAggregate(1.31, 40, 80, 2);
	toFile(*hex);
}

void testHexagonalAggregateRot(double b, double g)
{
	Particle *hex = new HexagonalAggregate(1.31, 40, 80, 2);
	double beta = (b*M_PI)/180;
	double gamma = (g*M_PI)/180;
	hex->Rotate(beta, gamma, 0);
	toFile(*hex);
}
//void testCompareParticles()
//{
//	Particle *hex1 = new TiltedHexagonal(40, 100, 1.31, 0);
//	Particle *hex2 = new Hexagonal(40, 100, 1.31);
//	int gg = 0;
//}
