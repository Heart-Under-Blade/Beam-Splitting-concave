#include "RegularIncidence.h"

#include "Beam.h"
#include "Splitting.h"

void RegularIncidence::ComputeDirections(Beam &beam, Splitting &splitter) const
{
	splitter.ComputePolarisationParams(beam);

	if (beam.isInside)
	{
		splitter.ComputeReflectedDirection(splitter.inBeam.direction);
		splitter.ComputeRefractedDirection(splitter.outBeam.direction);
	}
	else
	{
		splitter.ComputeReflectedDirection(splitter.outBeam.direction);
		splitter.ComputeRefractedDirection(splitter.inBeam.direction);
	}
#ifdef _DEBUG // DEB
	splitter.inBeam.dirs.push_back(splitter.inBeam.direction);
	splitter.outBeam.dirs.push_back(splitter.outBeam.direction);
#endif
	splitter.inBeam.polarizationBasis = beam.polarizationBasis;
	splitter.outBeam.polarizationBasis = beam.polarizationBasis;
}

void RegularIncidence::ComputeJonesMatrices(Beam &beam, Splitting &splitter) const
{
	splitter.inBeam.J = beam.J;
	splitter.outBeam.J = beam.J;

	if (beam.isInside)
	{
		double cosG = Point3f::DotProduct(splitter.m_normal, splitter.outBeam.direction);

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
		double cosB = Point3f::DotProduct(splitter.m_normal, splitter.inBeam.direction);

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
