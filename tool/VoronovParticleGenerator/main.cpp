#include <iostream>

#include "VoronoiLattice.h"

using namespace std;

int main()
{
	try {
//		VoronoiLattice(200, 3);
		VoronoiLattice(200, 5);
	} catch (...) {
		cout << "Fuck!";
	}
	cout << "Done.";
	return 0;
}
