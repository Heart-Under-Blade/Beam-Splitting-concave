#include "RegularIncidence.h"

#include "Beam.h"
#include "Splitting.h"

RegularIncidence::RegularIncidence(const complex &ri)
	: Incidence(ri)
{
}

void RegularIncidence::ComputeDirections(Beam &beam,
										 SplittedBeams<Beam> &beams) const
{
	beam.RotateJones(m_normal);

	if (beam.isInside)
	{
		ComputeReflectedDirection(beams.internal.direction);
		ComputeRefractedDirection(beams.external.direction);
	}
	else
	{
		ComputeReflectedDirection(beams.external.direction);
		ComputeRefractedDirection(beams.internal.direction);
	}
#ifdef _DEBUG // DEB
	beams.internal.dirs.push_back(beams.internal.direction);
	beams.external.dirs.push_back(beams.external.direction);
#endif
	beams.internal.polarizationBasis = beam.polarizationBasis;
	beams.external.polarizationBasis = beam.polarizationBasis;
}

void RegularIncidence::ComputeJonesMatrices(Beam &beam,
											SplittedBeams<Beam> &beams) const
{
	beams.internal.Jones = beam.Jones;
	beams.external.Jones = beam.Jones;

	if (beam.isInside)
	{
		double cosG = Point3f::DotProduct(normal, beams.external.direction);

		complex tmp0 = m_ri * cosA;
		complex tmp1 = m_ri * cosG;

		complex Tv0 = tmp1 + cosA;
		complex Th0 = tmp0 + cosG;

		complex tmp = 2.0 * tmp0;
		beams.external.MultiplyJonesMatrix(tmp/Tv0, tmp/Th0);

		complex Tv = cosA - tmp1;
		complex Th = tmp0 - cosG;
		beams.internal.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);
	}
	else
	{
		double cosB = Point3f::DotProduct(normal, beams.internal.direction);

		complex tmp0 = m_ri * cosA;
		complex tmp1 = m_ri * cosB;

		complex Tv0 = tmp0 + cosB;
		complex Th0 = tmp1 + cosA;

		complex Tv = tmp0 - cosB;
		complex Th = cosA - tmp1;
		beams.external.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);

		double cos2A = 2.0*cosA;
		beams.internal.MultiplyJonesMatrix(cos2A/Tv0, cos2A/Th0);
	}
}

void Incidence::ComputeOpticalPaths(const PathedBeam &parentBeam,
									SplittedBeams<PathedBeam> &beams) const
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
		double path = parentBeam.ComputeSegmentOpticalPath(reRiEff,
														   beams.internal.Center());
#ifdef _DEBUG // DEB
		beams.internal.ops = parentBeam.ops;
		beams.external.ops = parentBeam.ops;
#endif
		path += parentBeam.opticalPath;
		beams.internal.AddOpticalPath(path);
		beams.external.AddOpticalPath(path);
	}
}

