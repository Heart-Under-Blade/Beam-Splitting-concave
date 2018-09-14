#include "NormalIncidence.h"

#include "Beam.h"
#include "Splitting.h"

void NormalIncidence::ComputeLightParams(const Beam &incidentBeam,
										 Splitting &splitter) const
{
	const Point3f &dir = (incidentBeam.location == Location::In)
			? incidentBeam.direction
			: -incidentBeam.direction;

	splitter.inBeam.SetLight(-dir, incidentBeam.polarizationBasis);
	splitter.outBeam.SetLight(dir, incidentBeam.polarizationBasis);
}

void NormalIncidence::ComputeJonesMatrices(const Beam &incidentBeam,
										   Splitting &splitter) const
{
	splitter.inBeam.J = incidentBeam.J;
	splitter.outBeam.J = incidentBeam.J;

	complex f;

	if (incidentBeam.location == Location::In)
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
