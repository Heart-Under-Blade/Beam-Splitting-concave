#include "RegularIncidence.h"

#include "Beam.h"
#include "Splitting.h"

void RegularIncidence::ComputeDirections(Beam &beam, SplittedBeams<Beam> &beams,
										 Splitting &splitter) const
{
	beam.RotateJones(splitter.m_normal);

	if (beam.isInside)
	{
		splitter.ComputeReflectedDirection(beams.internal.direction);
		splitter.ComputeRefractedDirection(beams.external.direction);
	}
	else
	{
		splitter.ComputeReflectedDirection(beams.external.direction);
		splitter.ComputeRefractedDirection(beams.internal.direction);
	}

	beams.internal.polarizationBasis = beam.polarizationBasis;
	beams.external.polarizationBasis = beam.polarizationBasis;
}

void RegularIncidence::ComputeJonesMatrices(Beam &beam,
											SplittedBeams<Beam> &beams,
											Splitting &splitter) const
{
	beams.internal.Jones = beam.Jones;
	beams.external.Jones = beam.Jones;

	if (beam.isInside)
	{
		double cosG = Point3f::DotProduct(splitter.m_normal, beams.external.direction);

		complex tmp0 = splitter.m_ri * splitter.cosA;
		complex tmp1 = splitter.m_ri * cosG;

		complex Tv0 = tmp1 + splitter.cosA;
		complex Th0 = tmp0 + cosG;

		complex tmp = 2.0 * tmp0;
		beams.external.MultiplyJonesMatrix(tmp/Tv0, tmp/Th0);

		complex Tv = splitter.cosA - tmp1;
		complex Th = tmp0 - cosG;
		beams.internal.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);
	}
	else
	{
		double cosB = Point3f::DotProduct(splitter.m_normal, beams.internal.direction);

		complex tmp0 = splitter.m_ri * splitter.cosA;
		complex tmp1 = splitter.m_ri * cosB;

		complex Tv0 = tmp0 + cosB;
		complex Th0 = tmp1 + splitter.cosA;

		complex Tv = tmp0 - cosB;
		complex Th = splitter.cosA - tmp1;
		beams.external.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);

		double cos2A = 2.0*splitter.cosA;
		beams.internal.MultiplyJonesMatrix(cos2A/Tv0, cos2A/Th0);
	}
}

