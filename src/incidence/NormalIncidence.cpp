#include "NormalIncidence.h"

#include "Beam.h"
#include "Splitting.h"

NormalIncidence::NormalIncidence()
{
}

void NormalIncidence::SetSplitting(Splitting *splitting)
{
	Incidence::SetSplitting(splitting);
	complex tmp = (m_splitting->m_ri + 1.0);
	fresnels[0] = (2.0 * m_splitting->m_ri)/tmp;
	fresnels[1] = (1.0 - m_splitting->m_ri)/tmp;
	fresnels[2] = 2.0/tmp;
	fresnels[3] = (m_splitting->m_ri - 1.0)/tmp;
}

void NormalIncidence::ComputeDirections(Beam &beam, BeamPair<Beam> &beams,
										bool isBeamInside) const
{
	beams.internal.direction = (isBeamInside) ? -beam.direction
												: beam.direction;

	beams.external.direction = (isBeamInside) ? beam.direction
											   : -beam.direction;

	beams.internal.polarizationBasis = beam.polarizationBasis;
	beams.external.polarizationBasis = beam.polarizationBasis;
}

void NormalIncidence::ComputeJonesMatrices(Beam &beam, BeamPair<Beam> &beams,
										   bool isBeamInside) const
{
	beams.internal.Jones = beam.Jones;
	beams.external.Jones = beam.Jones;

	if (isBeamInside)
	{
		beams.external.MultiplyByFresnel(fresnels[0], fresnels[0]);
		beams.internal.MultiplyByFresnel(fresnels[1], -fresnels[1]);
	}
	else
	{
		beams.external.MultiplyByFresnel(fresnels[2], fresnels[2]);
		beams.internal.MultiplyByFresnel(fresnels[3], -fresnels[3]);
	}
}
