#include "RegularIncidence.h"

#include "Beam.h"
#include "Splitting.h"

void RegularIncidence::ComputeLightParams(const Beam &incidentBeam,
										  Splitting &splitter) const
{
	if (incidentBeam.location == Location::In)
	{
		Point3f reflDir = splitter.r - splitter.m_normal;
		Normalize(reflDir);
		splitter.inBeam.SetLight(reflDir, incidentBeam.polarizationBasis);

		Point3f refrDir = splitter.r/sqrt(splitter.s) + splitter.m_normal;
		Normalize(refrDir);
		splitter.outBeam.SetLight(refrDir, incidentBeam.polarizationBasis);
	}
	else
	{
		splitter.inBeam.polarizationBasis = incidentBeam.polarizationBasis;

		Point3f refrDir;
		const Point3f &dir = incidentBeam.direction;

		Point3f r = dir/splitter.cosA + splitter.m_normal;

		refrDir = r + splitter.m_normal;
		Normalize(refrDir);

		splitter.ComputeInternalRefractiveDirection(splitter.r, splitter.m_normal,
													splitter.inBeam.direction);

		splitter.outBeam.SetLight(refrDir, incidentBeam.polarizationBasis);

	}
}

void RegularIncidence::ComputeJonesMatrices(const Beam &incidentBeam,
											Splitting &splitter) const

{
	splitter.inBeam.J = incidentBeam.J;
	splitter.outBeam.J = incidentBeam.J;

	if (incidentBeam.location == Location::In)
	{
		double cosG = DotProduct(splitter.m_normal, splitter.outBeam.direction);

		complex tmp0 = splitter.m_ri * splitter.cosA;
		complex tmp1 = splitter.m_ri * cosG;

		complex Tv0 = tmp1 + splitter.cosA;
		complex Th0 = tmp0 + cosG;

		complex tmp = 2.0 * tmp0;
		splitter.outBeam.MultiplyJonesMatrix(tmp/Tv0, tmp/Th0);

		complex Tv = splitter.cosA - tmp1;
		complex Th = tmp0 - cosG;
		splitter.inBeam.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);
	}
	else
	{
		double cosB = DotProduct(-splitter.m_normal, splitter.inBeam.direction);

		complex tmp0 = splitter.m_ri * splitter.cosA;
		complex tmp1 = splitter.m_ri * cosB;

		complex Tv0 = tmp0 + cosB;
		complex Th0 = tmp1 + splitter.cosA;

		complex Tv = tmp0 - cosB;
		complex Th = splitter.cosA - tmp1;
		splitter.outBeam.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);

		double cos2A = 2.0*splitter.cosA;
		splitter.inBeam.MultiplyJonesMatrix(cos2A/Tv0, cos2A/Th0);
	}
}
