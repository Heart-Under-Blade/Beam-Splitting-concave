#include "ScatteringConvex.h"

ScatteringConvex::ScatteringConvex(Particle *particle, const Light &incidentLight,
								   int maxActNo)
	: Scattering(particle, incidentLight, maxActNo)
{
}

void ScatteringConvex::PushBeamsToBuffer(Facet *facet, BeamPair<Beam> &beams,
										 std::vector<Beam> &scatteredBeams)
{
<<<<<<< HEAD
//	m_particle->Rotate(beta, gamma, 0);
=======
	Track tr = m_originalBeam;
	tr.Update(facet);
	tr.RecomputeTrackId(0, facet->index);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

	beams.external.CopyTrack(tr);
	beams.external.SetLocation(false);
	scatteredBeams.push_back(beams.external);

	beams.internal.CopyTrack(tr);
	beams.internal.SetLocation(true);

	if (IsTerminalAct(beams.internal))
	{
		if (!beams.internal.isInside)
		{
			ReleaseBeam(beams.internal);
		}
	}
	else
	{
		PushBeamToTree(beams.internal);
	}

#ifdef _CHECK_ENERGY_BALANCE
	ComputeFacetEnergy(facet->in_normal, beams.external);
#endif
}

void ScatteringConvex::SplitOriginalBeam(std::vector<Beam> &externalBeams)
{
	m_visibleFacets.nElems = 0;
	FindVisibleFacets(m_originalBeam, m_lightChecker, 0, m_particle->nElems, m_visibleFacets);

	for (int i = 0; i < m_visibleFacets.nElems; ++i)
	{
		Facet *facet = m_visibleFacets.elems[i];

		m_splitting.SetBeams(*facet);
		ComputeOpticalBeamParams(facet, m_originalBeam);
		PushBeamsToBuffer(facet, m_splitting.beams, externalBeams);
	}
}

void ScatteringConvex::SelectVisibleFacets(const Beam &beam, Array<Facet*> &facets)
{
<<<<<<< HEAD
	Beam outBeam;
	const Point3f &incidentDir = incidentBeam.direction;

	// ext. normal uses in this calculating
	const Point3f &normal = m_facets[facetID].ex_normal;
	m_splitting.ComputeCosA(normal, incidentDir);

	if (!m_splitting.IsIncident())
	{
		return false;
	}

#ifdef _DEBUG // DEB
	if (outBeams.size() == 329)
		int ddd = 0;
#endif
	Intersect(facetID, incidentBeam, outBeam);

	if (outBeam.nVertices < MIN_VERTEX_NUM)
	{
		return false;
	}

	inBeam = outBeam;
	m_splitting.ComputeSplittingParams(incidentBeam.direction, normal);

	if (!m_splitting.IsNormalIncidence())
	{	// regular incidence
		ComputePolarisationParams(incidentBeam.direction, normal, incidentBeam);

		if (!m_splitting.IsCompleteReflection())
		{
			outBeam.id = RecomputeTrackId(incidentBeam.id, facetID);

			m_splitting.ComputeRegularBeamsParams(normal, incidentBeam,
												  inBeam, outBeam);
			outBeam.nActs = incidentBeam.nActs + 1;
			outBeam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
			outBeam.lastFacetId = facetID;
			outBeams.push_back(outBeam);
#ifdef _DEBUG // DEB
			if (outBeams.size() == 330)
				int ddd = 0;
#endif
		}
		else // complete internal reflection incidence
		{
			m_splitting.ComputeCRBeamParams(normal, incidentBeam, inBeam);
		}
	}
	else
	{	// normal incidence
		m_splitting.ComputeNormalBeamParams(incidentBeam, inBeam, outBeam);

		outBeam.nActs = incidentBeam.nActs + 1;
		outBeam.id = RecomputeTrackId(incidentBeam.id, facetID);
		double path = m_splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
		outBeam.opticalPath += path;
		outBeam.lastFacetId = facetID;
		outBeams.push_back(outBeam);
#ifdef _DEBUG // DEB
			if (outBeams.size() == 330)
				int ddd = 0;
#endif
	}

	return true;
=======
	FindVisibleFacets(beam, m_lightChecker, 0, m_particle->nElems, facets);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
}
