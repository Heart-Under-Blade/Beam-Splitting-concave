#include "TotalReflectionIncidence.h"

#include "Beam.h"
#include "Splitting.h"

TotalReflectionIncidence::TotalReflectionIncidence()
{
}

void TotalReflectionIncidence::ComputeDirections(Beam &beam,
												 BeamPair<Beam> &beams,
												 bool isBeamInside) const
{
	m_splitting->ComputeReflectedDirection(beams.internal.direction);
	beams.internal.polarizationBasis = beam.polarizationBasis;
}

void TotalReflectionIncidence::ComputeJonesMatrices(Beam &parentBeam,
													BeamPair<Beam> &beams,
													bool isBeamInside) const
{
	const double bf = m_splitting->reRiEff*(1.0 - m_splitting->cosA2) - 1.0;
	double im = (bf > 0) ? sqrt(bf) : 0;

	const complex sq(0, im);
	complex tmp0 = m_splitting->m_ri * m_splitting->cosA;
	complex tmp1 = m_splitting->m_ri * sq;

	complex cv = (m_splitting->cosA - tmp1)/(tmp1 + m_splitting->cosA);
	complex ch = (tmp0 - sq)/(tmp0 + sq);

	beams.internal.Jones = parentBeam.Jones;
	beams.internal.MultiplyByFresnel(cv, ch);
}
