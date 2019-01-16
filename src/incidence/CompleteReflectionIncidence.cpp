#include "CompleteReflectionIncidence.h"

#include "Beam.h"
#include "Splitting.h"

CompleteReflectionIncidence::CompleteReflectionIncidence(const complex &ri)
	: Incidence(ri)
{
}

void CompleteReflectionIncidence::ComputeDirections(Beam &beam,
													SplittedBeams<Beam> &beams) const
{
	ComputeReflectedDirection(beams.internal.direction);
	beams.internal.polarizationBasis = beam.polarizationBasis;
#ifdef _DEBUG // DEB
	beams.internal.dirs.push_back(beams.internal.direction);
#endif
}

void CompleteReflectionIncidence::ComputeJonesMatrices(Beam &parentBeam,
													   SplittedBeams<Beam> &beams) const
{
	const double bf = reRiEff*(1.0 - cosA2) - 1.0;
    double im = (bf > 0) ? sqrt(bf) : 0;

    const complex sq(0, im);
	complex tmp0 = m_ri * cosA;
	complex tmp1 = m_ri * sq;

	complex cv = (cosA - tmp1)/(tmp1 + cosA);
    complex ch = (tmp0 - sq)/(tmp0 + sq);

	beams.internal.Jones = parentBeam.Jones;
	beams.internal.MultiplyJonesMatrix(cv, ch);
}

void CompleteReflectionIncidence::ComputeOpticalPaths(const PathedBeam &incidentBeam,
													  SplittedBeams<PathedBeam> &beams) const
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
		double path = incidentBeam.ComputeSegmentOpticalPath(reRiEff,
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
