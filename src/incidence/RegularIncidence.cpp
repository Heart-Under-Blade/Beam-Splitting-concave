#include "RegularIncidence.h"

#include "Beam.h"
#include "Splitting.h"

RegularIncidence::RegularIncidence()
{
}

void RegularIncidence::ComputeDirections(Beam &beam,
										 BeamPair<Beam> &beams) const
{
	beam.RotateJones(m_splitting->facetNormal);

	if (beam.isInside)
	{
		m_splitting->ComputeReflectedDirection(beams.internal.direction);
		m_splitting->ComputeRefractedDirection(beams.external.direction);
	}
	else
	{
		m_splitting->ComputeReflectedDirection(beams.external.direction);
		m_splitting->ComputeRefractedDirection(beams.internal.direction);
	}

	beams.internal.polarizationBasis = beam.polarizationBasis;
	beams.external.polarizationBasis = beam.polarizationBasis;
}

void RegularIncidence::ComputeJonesMatrices(Beam &beam,
											BeamPair<Beam> &beams) const
{
	beams.internal.Jones = beam.Jones;
	beams.external.Jones = beam.Jones;

	if (beam.isInside)
	{
		double cosG = Point3f::DotProduct(m_splitting->facetNormal, beams.external.direction);

		complex tmp0 = m_splitting->m_ri * m_splitting->cosA;
		complex tmp1 = m_splitting->m_ri * cosG;

		complex Tv0 = tmp1 + m_splitting->cosA;
		complex Th0 = tmp0 + cosG;

		complex tmp = 2.0 * tmp0;
		beams.external.MultiplyJonesMatrix(tmp/Tv0, tmp/Th0);

		complex Tv = m_splitting->cosA - tmp1;
		complex Th = tmp0 - cosG;
		beams.internal.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);
	}
	else
	{
		double cosB = Point3f::DotProduct(m_splitting->facetNormal, beams.internal.direction);

		complex tmp0 = m_splitting->m_ri * m_splitting->cosA;
		complex tmp1 = m_splitting->m_ri * cosB;

		complex Tv0 = tmp0 + cosB;
		complex Th0 = tmp1 + m_splitting->cosA;

		complex Tv = tmp0 - cosB;
		complex Th = m_splitting->cosA - tmp1;
		beams.external.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);

		double cos2A = 2.0*m_splitting->cosA;
		beams.internal.MultiplyJonesMatrix(cos2A/Tv0, cos2A/Th0);
	}
}
