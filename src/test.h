/**
  File for various tests
*/

#pragma once

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

void testHexagonBuilding()
{
	Particle hex = Hexagonal(20, 100, 1.31);
	outputParticle(hex);
}

void testConcaveHexagonBuilding()
{
	Particle *hex = new ConcaveHexagonal(20, 100, 1.31, 10);
	hex->Rotate(1.5708/*90 degrees*/, 0, 0);
	outputParticle(*hex);
}

void testHexagonRotate()
{
	Particle *hex = new Hexagonal(20, 100, 1.31);
	hex->Rotate(1.5708/*90 degrees*/, 0, 0);

	outputParticle(*hex);
}


