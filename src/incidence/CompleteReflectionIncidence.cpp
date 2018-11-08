#include "CompleteReflectionIncidence.h"

#include "Beam.h"
#include "Splitting.h"

void CompleteReflectionIncidence::ComputeDirections(Beam &beam,
													SplittedBeams<Beam> &beams,
													Splitting &splitter) const
{
	beam.RotateJones(splitter.m_normal);
	splitter.ComputeReflectedDirection(beams.internal.direction);
	beams.internal.polarizationBasis = beam.polarizationBasis;
}

void CompleteReflectionIncidence::ComputeJonesMatrices(Beam &parentBeam,
													   SplittedBeams<Beam> &beams,
													   Splitting &splitter) const
{
    const double bf = splitter.reRiEff*(1.0 - splitter.cosA*splitter.cosA) - 1.0;
    double im = (bf > 0) ? sqrt(bf) : 0;

    const complex sq(0, im);
    complex tmp0 = splitter.m_ri * splitter.cosA;
    complex tmp1 = splitter.m_ri * sq;

    complex cv = (splitter.cosA - tmp1)/(tmp1 + splitter.cosA);
    complex ch = (tmp0 - sq)/(tmp0 + sq);

	beams.internal.Jones = parentBeam.Jones;
	beams.internal.MultiplyJonesMatrix(cv, ch);
}

void CompleteReflectionIncidence::ComputeOpticalPaths(const PathedBeam &incidentBeam,
													  SplittedBeams<PathedBeam> &beams,
													  Splitting &splitter) const
{
	if (incidentBeam.opticalPath < FLT_EPSILON)
	{
		Point3f p = beams.internal.Center();
		double path = incidentBeam.ComputeIncidentOpticalPath(incidentBeam.direction, p);
		beams.internal.opticalPath = 0;
		beams.internal.AddOpticalPath(path);
	}
	else
	{
		double path = incidentBeam.ComputeSegmentOpticalPath(splitter.reRiEff,
															 beams.internal.Center());
		path += incidentBeam.opticalPath;
#ifdef _DEBUG // DEB
		beams.internal.ops = incidentBeam.ops;
		beams.internal.ops.push_back(path);
#endif
		path += incidentBeam.opticalPath;
		beams.internal.AddOpticalPath(path);
	}
}
