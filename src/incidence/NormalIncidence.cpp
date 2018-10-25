#include "NormalIncidence.h"

#include "Beam.h"
#include "Splitting.h"

void NormalIncidence::ComputeDirections(Beam &beam, Splitting &splitter) const
{
	splitter.inBeam.direction = (beam.isInside) ? -beam.direction
												: beam.direction;

	splitter.outBeam.direction = (beam.isInside) ? beam.direction
												 : -beam.direction;

	splitter.inBeam.polarizationBasis = beam.polarizationBasis;
	splitter.outBeam.polarizationBasis = beam.polarizationBasis;

#ifdef _DEBUG // DEB
	splitter.inBeam.dirs.push_back(splitter.inBeam.direction);
	splitter.outBeam.dirs.push_back(splitter.outBeam.direction);
#endif
}

void NormalIncidence::ComputeJonesMatrices(Beam &beam, Splitting &splitter) const
{
	splitter.inBeam.Jones = beam.Jones;
	splitter.outBeam.Jones = beam.Jones;

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
