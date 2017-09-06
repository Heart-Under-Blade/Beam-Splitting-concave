#include "Hexagonal.h"
#include "ConcaveHexagonal.h"
#include "HexagonalAggregate.h"
#include "BulletRosette.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>

void writeFacet(int beg, int end, std::ofstream &file, Particle *particle)
{
	for (int i = beg; i < end; ++i)
	{
		for (int j = particle->facets[i].size-1; j >= 0; --j)
		{
			file << particle->facets[i].arr[j].cx << ' '
				 << particle->facets[i].arr[j].cy << ' '
				 << particle->facets[i].arr[j].cz << ' '
				 << '\n';
		}
	}
}

int main()
{
	int type;
	double height, diameter;

	printf("Particle type \n0 - Hex prism;\n1 - Concave hex prism;\n2 - Hex prism aggregate\n"\
		   "3 - Bullet-rosette\n");
	scanf("%d", &type);

	printf("\nHeight:\n");
	scanf("%lf", &height);

	printf("\nDiameter:\n");
	scanf("%lf", &diameter);

	Particle *particle;

	if (type == 1)
	{
		double cavityDept;
		printf("\nCavity dept:\n");
		scanf("%lf", &cavityDept);

		particle = new ConcaveHexagonal(1.31, diameter, height, cavityDept);
	}
	else if (type == 0)
	{
		particle = new Hexagonal(1.31, diameter, height);
	}
	else if (type == 2)
	{
		particle = new HexagonalAggregate(1.31, diameter, height, 2);
	}
	else if (type == 3)
	{
		double sup = (diameter*sqrt(3)*1,8807264653463320123608375958293)/4;
		particle = new BulletRosette(1.31, diameter, height, sup);
	}

	// output file
	std::ofstream file("hexagon_crystal.crystal", std::ios::out);
	file << particle->facetNum << "\n";

	if (type == 3)
	{
		for (int i = 0; i < particle->facetNum; ++i)
		{
			file << particle->facets[i].size << "\n";
		}
	}
	else
	{
		std::vector<int> vNums;
		for (int i = 0; i < particle->facetNum; ++i)
		{
			vNums.push_back(particle->facets[i].size);
		}

		std::sort(vNums.begin(), vNums.end(), std::greater<int>());

		for (auto n : vNums)
		{
			file << n << "\n";
		}
	}

	if (type == 1)
	{
		writeFacet(6, 12, file, particle);
		writeFacet(0, 6, file, particle);
		writeFacet(12, 18, file, particle);
	}
	else if (type == 0)
	{
		writeFacet(0, 1, file, particle);
		writeFacet(7, 8, file, particle);
		writeFacet(1, 7, file, particle);
	}
	else if (type == 2)
	{
		writeFacet(0, 1, file, particle);
		writeFacet(7, 8, file, particle);
		writeFacet(8, 9, file, particle);
		writeFacet(15, 16, file, particle);
		writeFacet(1, 7, file, particle);
		writeFacet(9, 15, file, particle);
	}
	else if (type == 3)
	{
		writeFacet(0, 78, file, particle);
	}

	printf("\nAll done. Press anykey...");
	getchar();
	return 0;
}
