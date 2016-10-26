#include "Hexagonal.h"
#include "ConcaveHexagonal.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>

void writeFacet(int beg, int end, std::ofstream &file, Particle *particle)
{
	for (int i = beg; i < end; ++i)
	{
		for (int j = particle->vertexNums[i]-1; j >= 0; --j)
		{
			file << particle->facets[i][j].cx << ' '
				 << particle->facets[i][j].cy << ' '
				 << particle->facets[i][j].cz << ' '
				 << '\n';
		}
	}
}

int main()
{
	int type;
	double halfHeight, radius;

	printf("Particle type (0 - Prism, 1 - Concave prism):\n");
	scanf("%d", &type);

	printf("\nHalfheight:\n");
	scanf("%lf", &halfHeight);

	printf("\nRadius:\n");
	scanf("%lf", &radius);

	Particle *particle;

	if (type == 1)
	{
		double cavityDept;
		printf("\nCavity dept:\n");
		scanf("%lf", &cavityDept);

		particle = new ConcaveHexagonal(radius, halfHeight, 1.31, cavityDept);
	}
	else
	{
		particle = new Hexagonal(radius, halfHeight, 1.31);
	}

	// output file
	std::ofstream file("hexagon_crystal.crystal", std::ios::out);
	file << particle->facetNum << "\n";

	std::vector<int> vNums;
	for (int i = 0; i < particle->facetNum; ++i)
	{
		vNums.push_back(particle->vertexNums[i]);
	}

	std::sort(vNums.begin(), vNums.end(), std::greater<int>());

	for (auto n : vNums)
	{
		file << n << "\n";
	}

	if (type == 1)
	{
		writeFacet(6, 12, file, particle);
		writeFacet(0, 6, file, particle);
		writeFacet(12, 18, file, particle);
	}
	else
	{
		writeFacet(0, 1, file, particle);
		writeFacet(7, 8, file, particle);
		writeFacet(1, 7, file, particle);
	}

	printf("\nAll done. Press anykey...");
	getchar();
	return 0;
}
