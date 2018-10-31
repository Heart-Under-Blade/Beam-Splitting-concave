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
#ifdef _DEBUG // DEB
	beams.internal.dirs.push_back(beams.internal.direction);
	beams.external.dirs.push_back(beams.external.direction);
#endif
	beams.internal.polarizationBasis = beam.polarizationBasis;
	beams.external.polarizationBasis = beam.polarizationBasis;
}

void RegularIncidence::ComputeJonesMatrices(Beam &beam, Splitting &splitter)
{
	splitter.internal.Jones = beam.Jones;
	splitter.external.Jones = beam.Jones;

	if (beam.isInside)
	{
		double cosG = Point3f::DotProduct(splitter.m_normal, splitter.external.direction);

		complex tmp0 = splitter.m_ri * splitter.cosA;
		complex tmp1 = splitter.m_ri * cosG;

		complex Tv0 = tmp1 + splitter.cosA;
		complex Th0 = tmp0 + cosG;

		complex tmp = 2.0 * tmp0;
		splitter.external.MultiplyJonesMatrix(tmp/Tv0, tmp/Th0);

		complex Tv = splitter.cosA - tmp1;
		complex Th = tmp0 - cosG;
		splitter.internal.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);
	}
	else
	{
		double cosB = Point3f::DotProduct(splitter.m_normal, splitter.internal.direction);

		complex tmp0 = splitter.m_ri * splitter.cosA;
		complex tmp1 = splitter.m_ri * cosB;

		complex Tv0 = tmp0 + cosB;
		complex Th0 = tmp1 + splitter.cosA;

		complex Tv = tmp0 - cosB;
		complex Th = splitter.cosA - tmp1;
		splitter.external.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);

		double cos2A = 2.0*splitter.cosA;
		splitter.internal.MultiplyJonesMatrix(cos2A/Tv0, cos2A/Th0);
	}
}


void Incidence::ComputeOpticalPaths(const PathedBeam &parentBeam,
									SplittedBeams<PathedBeam> &beams,
									Splitting &splitter) const
{
	if (parentBeam.opticalPath < FLT_EPSILON)
	{
		Point3f p = beams.internal.Center();
		double path = parentBeam.ComputeIncidentOpticalPath(parentBeam.direction, p);
		beams.internal.AddOpticalPath(path);
		beams.external.AddOpticalPath(path);
	}
	else
	{
		double path = parentBeam.ComputeSegmentOpticalPath(splitter.reRiEff,
													 beams.internal.Center());
#ifdef _DEBUG // DEB
		splitter.internal.ops = beam.ops;
		splitter.outBeam.ops = beam.ops;
#endif
		path += parentBeam.opticalPath;
		beams.internal.AddOpticalPath(path);
		beams.external.AddOpticalPath(path);
	}
}

