#include "NormalIncidence.h"

#include "Beam.h"
#include "Splitting.h"

void NormalIncidence::ComputeDirections(Beam &beam, SplittedBeams<Beam> &beams,
										Splitting &splitter) const
{
	beams.internal.direction = (beam.isInside) ? -beam.direction
												: beam.direction;

	beams.external.direction = (beam.isInside) ? beam.direction
												 : -beam.direction;

	beams.internal.polarizationBasis = beam.polarizationBasis;
	beams.external.polarizationBasis = beam.polarizationBasis;
}

void NormalIncidence::ComputeJonesMatrices(Beam &beam, SplittedBeams<Beam> &beams,
										   Splitting &splitter) const
{
	beams.internal.Jones = beam.Jones;
	beams.external.Jones = beam.Jones;

	complex f;

	if (beam.isInside)
	{
		f = (2.0 * splitter.m_ri)/(1.0 + splitter.m_ri); // OPT: вынести целиком
		beams.external.MultiplyJonesMatrix(f, f);

		f = (1.0 - splitter.m_ri)/(1.0 + splitter.m_ri); // OPT: вынести целиком
		beams.internal.MultiplyJonesMatrix(f, -f);
	}
	else
	{
		f = 2.0/(splitter.m_ri + 1.0); // OPT: вынести целиком
		beams.external.MultiplyJonesMatrix(f, f);

		f = (splitter.m_ri - 1.0)/(splitter.m_ri + 1.0); // OPT: вынести целиком
		beams.internal.MultiplyJonesMatrix(f, -f);
	}
}
