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

void NormalIncidence::ComputeDirections(Beam &beam, BeamPair<Beam> &beams)
{
	beams.internal.direction = (beam.isInside) ? -beam.direction
												: beam.direction;

	beams.external.direction = (beam.isInside) ? beam.direction
											   : -beam.direction;

	beams.internal.polarizationBasis = beam.polarizationBasis;
	beams.external.polarizationBasis = beam.polarizationBasis;

#ifdef _DEBUG // DEB
	beams.internal.dirs.push_back(beams.internal.direction);
	beams.external.dirs.push_back(beams.external.direction);
#endif
}

void NormalIncidence::ComputeJonesMatrices(Beam &beam, BeamPair<Beam> &beams)
{
	beams.internal.Jones = beam.Jones;
	beams.external.Jones = beam.Jones;

	if (beam.isInside)
	{
		beams.external.MultiplyJonesMatrix(fresnels[0], fresnels[0]);
		beams.internal.MultiplyJonesMatrix(fresnels[1], -fresnels[1]);
	}
	else
	{
		beams.external.MultiplyJonesMatrix(fresnels[2], fresnels[2]);
		beams.internal.MultiplyJonesMatrix(fresnels[3], -fresnels[3]);
	}
}
