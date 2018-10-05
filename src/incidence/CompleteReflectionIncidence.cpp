#include "CompleteReflectionIncidence.h"

#include "Beam.h"
#include "Splitting.h"

void CompleteReflectionIncidence::ComputeDirections(const Beam &incidentBeam,
                                                     Splitting &splitter) const
{
	splitter.ComputeReflectedDirection(splitter.inBeam.direction);
	splitter.inBeam.polarizationBasis = incidentBeam.polarizationBasis;
}

void CompleteReflectionIncidence::ComputeJonesMatrices(const Beam &incidentBeam,
                                                       Splitting &splitter) const
{
    const double bf = splitter.reRiEff*(1.0 - splitter.cosA*splitter.cosA) - 1.0;
    double im = (bf > 0) ? sqrt(bf) : 0;

    const complex sq(0, im);
    complex tmp0 = splitter.m_ri * splitter.cosA;
    complex tmp1 = splitter.m_ri * sq;

    complex cv = (splitter.cosA - tmp1)/(tmp1 + splitter.cosA);
    complex ch = (tmp0 - sq)/(tmp0 + sq);

    splitter.inBeam.J = incidentBeam.J;
    splitter.inBeam.MultiplyJonesMatrix(cv, ch);
}

void CompleteReflectionIncidence::ComputeOpticalPaths(const Beam &incidentBeam,
                                                      Splitting &splitter) const
{
	if (incidentBeam.opticalPath < FLT_EPSILON)
	{
		Point3f p = splitter.inBeam.Center();
		double path = splitter.ComputeIncidentOpticalPath(incidentBeam.direction, p);
		splitter.inBeam.opticalPath = 0;
		splitter.inBeam.AddOpticalPath(path);
	}
	else
	{
		double path = splitter.ComputeSegmentOpticalPath(incidentBeam,
														 splitter.inBeam.Center());
		path += incidentBeam.opticalPath;
#ifdef _DEBUG // DEB
		splitter.inBeam.ops = incidentBeam.ops;
		splitter.inBeam.ops.push_back(path);
#endif
		path += incidentBeam.opticalPath;
		splitter.inBeam.AddOpticalPath(path);
	}
}
