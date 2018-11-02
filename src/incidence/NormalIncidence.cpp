#include "NormalIncidence.h"

#include "Beam.h"
#include "Splitting.h"

void NormalIncidence::ComputeDirections(Beam &beam, Splitting &splitter)
{
	splitter.internal.direction = (beam.isInside) ? -beam.direction
												: beam.direction;

	splitter.external.direction = (beam.isInside) ? beam.direction
												 : -beam.direction;

	splitter.internal.polarizationBasis = beam.polarizationBasis;
	splitter.external.polarizationBasis = beam.polarizationBasis;

#ifdef _DEBUG // DEB
	splitter.internal.dirs.push_back(splitter.internal.direction);
	splitter.outBeam.dirs.push_back(splitter.outBeam.direction);
#endif
}

void NormalIncidence::ComputeJonesMatrices(Beam &beam, Splitting &splitter)
{
	splitter.internal.Jones = beam.Jones;
	splitter.external.Jones = beam.Jones;

	complex f;

	if (beam.isInside)
	{
		f = (2.0 * splitter.m_ri)/(1.0 + splitter.m_ri); // OPT: вынести целиком
		splitter.external.MultiplyJonesMatrix(f, f);

		f = (1.0 - splitter.m_ri)/(1.0 + splitter.m_ri); // OPT: вынести целиком
		splitter.internal.MultiplyJonesMatrix(f, -f);
	}
	else
	{
		f = 2.0/(splitter.m_ri + 1.0); // OPT: вынести целиком
		splitter.external.MultiplyJonesMatrix(f, f);

		f = (splitter.m_ri - 1.0)/(splitter.m_ri + 1.0); // OPT: вынести целиком
		splitter.internal.MultiplyJonesMatrix(f, -f);
	}
}
