#include "NormalIncidence.h"

#include "Beam.h"
#include "Splitting.h"

void NormalIncidence::ComputeDirections(Beam &beam, Splitting &splitter) const
{
	const Point3f &dir = (beam.isInside) ? beam.direction
										 : -beam.direction;

	splitter.inBeam.SetLight(-dir, beam.polarizationBasis);
	splitter.outBeam.SetLight(dir, beam.polarizationBasis);

#ifdef _DEBUG // DEB
	splitter.inBeam.dirs.push_back(splitter.inBeam.direction);
	splitter.outBeam.dirs.push_back(splitter.outBeam.direction);
#endif
}

void NormalIncidence::ComputeJonesMatrices(Beam &beam, Splitting &splitter) const
{
	splitter.inBeam.J = beam.J;
	splitter.outBeam.J = beam.J;

	complex f;

	if (beam.isInside)
	{
		f = (2.0 * splitter.m_ri)/(1.0 + splitter.m_ri); // OPT: вынести целиком
		splitter.outBeam.MultiplyJonesMatrix(f, f);

		f = (1.0 - splitter.m_ri)/(1.0 + splitter.m_ri); // OPT: вынести целиком
		splitter.inBeam.MultiplyJonesMatrix(f, -f);
	}
	else
	{
		f = 2.0/(splitter.m_ri + 1.0); // OPT: вынести целиком
		splitter.outBeam.MultiplyJonesMatrix(f, f);

		f = (splitter.m_ri - 1.0)/(splitter.m_ri + 1.0); // OPT: вынести целиком
		splitter.inBeam.MultiplyJonesMatrix(f, -f);
	}
}
