/**
  File for various tests
*/

#pragma once

#include <fstream>
#include <iostream>

#include "Hexagonal.h"
#include "ConcaveHexagonal.h"

void outputParticle(const Particle &particle)
{
	for (int i = 0; i < particle.facetNum; ++i)
	{
		std::cout << i << ": ";
		for (int j = 0; j < particle.vertexNums[i]; ++j)
		{
			Point3f p = particle.facets[i][j];
			std::cout << "("
					  << p.point[0] << ", "
					  << p.point[1] << ", "
					  << p.point[2]
					  << "), ";
		}

		std::cout << std::endl;
	}

	std::cout << std::endl << "Normals" << std::endl << std::endl;

	for (int i = 0; i < particle.facetNum; ++i)
	{
		std::cout << i << ": ";
		std::cout << "("
				  << particle.normals[i].point[0] << ", "
				  << particle.normals[i].point[1] << ", "
				  << particle.normals[i].point[2]
				  << "), ";
		std::cout << std::endl;
	}
}

void toFile(const Particle &particle)
{
	std::ofstream M("particle.dat", std::ios::out);

	for (int i = 0; i < particle.facetNum; ++i)
	{
		for (int j = 0; j < particle.vertexNums[i]; ++j)
		{
			Point3f p = particle.facets[i][j];
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
	Particle hex = Hexagonal(20, 100, 1.31);
	toFile(hex);
	outputParticle(hex);
}

void testConcaveHexagonRot()
{
	Particle *hex = new ConcaveHexagonal(20, 100, 1.31, 10);
	double beta = ((19 + 0.5)*M_PI)/(2.0*100);
	double gamma = ((36 + 0.5)*M_PI)/(3.0*101);
	hex->Rotate(beta, gamma, 0);
	toFile(*hex);
	outputParticle(*hex);
}

void testHexagonRotate()
{
	Particle *hex = new Hexagonal(20, 100, 1.31);
//	double beta = 19*M_PI/180;
//	double gamma = 90*M_PI/180;
	double beta = ((19 + 0.5)*M_PI)/(2.0*100);
	double gamma = ((36 + 0.5)*M_PI)/(3.0*101);
	hex->Rotate(beta, gamma, 0);

	toFile(*hex);
	outputParticle(*hex);
}


