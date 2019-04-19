#include "CompleteReflectionIncidence.h"

#include "Beam.h"
#include "Splitting.h"

CompleteReflectionIncidence::CompleteReflectionIncidence()
{
}

void CompleteReflectionIncidence::ComputeDirections(Beam &beam,
													BeamPair<Beam> &beams) const
{
	m_splitting->ComputeReflectedDirection(beams.internal.direction);
	beams.internal.polarizationBasis = beam.polarizationBasis;
}

void CompleteReflectionIncidence::ComputeJonesMatrices(Beam &parentBeam,
													   BeamPair<Beam> &beams) const
{
	const double bf = m_splitting->reRiEff*(1.0 - m_splitting->cosA2) - 1.0;
    double im = (bf > 0) ? sqrt(bf) : 0;

    const complex sq(0, im);
	complex tmp0 = m_splitting->m_ri * m_splitting->cosA;
	complex tmp1 = m_splitting->m_ri * sq;

	complex cv = (m_splitting->cosA - tmp1)/(tmp1 + m_splitting->cosA);
    complex ch = (tmp0 - sq)/(tmp0 + sq);

	beams.internal.Jones = parentBeam.Jones;
	beams.internal.MultiplyJonesMatrix(cv, ch);
}

void CompleteReflectionIncidence::ComputeOpticalPaths(const PathedBeam &incidentBeam,
													  BeamPair<PathedBeam> &beams) const
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
		double path = incidentBeam.ComputeSegmentOpticalPath(m_splitting->reRiEff,
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
