#include "ScatteringConvex.h"

<<<<<<< HEAD
#include <iostream>

ScatteringConvex::ScatteringConvex(Particle *particle, Light *incidentLight,
								   bool isOpticalPath, int nActs)
	: Scattering(particle, incidentLight, isOpticalPath, nActs)
{
}

void ScatteringConvex::ScatterLight(std::vector<Beam> &outBeams)
{
	m_incidentEnergy = 0;
	m_treeSize = 0;

	/// first extermal beam
	for (int facetID = 0; facetID < m_particle->nFacets; ++facetID)
=======
ScatteringConvex::ScatteringConvex(Particle *particle, const Light &incidentLight,
								   int maxActNo, const complex &refractiveIndex)
	: Scattering(particle, incidentLight, maxActNo, refractiveIndex)
{
}

void ScatteringConvex::SetWorkFacetsFromTracks()
{
	m_workFacets.nElems = 0;

	for (auto n : m_trackTreeNode->children)
	{
		auto facet = m_particle->GetActualFacet(n->m_facet->index);
		m_workFacets.Add(facet);
	}
}

void ScatteringConvex::ScatterLight(TrackNode *trackTree,
									std::vector<Beam> &scatteredBeams)
{
	m_hasTracks = true;
	m_trackTreeNode = trackTree;
	SetWorkFacetsFromTracks();

	SplitOriginalBeam(scatteredBeams);
	// TODO: to complete
}

void ScatteringConvex::PushBeamsToBuffer(Facet *facet, BeamPair<Beam> &beams,
										 std::vector<Beam> &scatteredBeams)
{
	Track tr = m_originalBeam;
	tr.Update(facet);
	tr.RecomputeTrackId(0, facet->index);

	beams.external.CopyTrack(tr);
	beams.external.SetLocation(false);
	scatteredBeams.push_back(beams.external);

	beams.internal.CopyTrack(tr);
	beams.internal.SetLocation(true);

	if (IsFinalAct(beams.internal))
>>>>>>> origin/refactor
	{
		const Point3f &inNormal = m_facets[facetID].in_normal;
		m_splitting.ComputeCosA(m_incidentDir, inNormal);

		if (!m_splitting.IsIncident()) /// beam is not incident to this facet
		{
			continue;
		}

		Beam inBeam, outBeam;
		SplitLightToBeams(facetID, inBeam, outBeam);

		auto newId = RecomputeTrackId(0, facetID);

		outBeam.id = newId;
		outBeam.lastFacetId = facetID;
		outBeam.nActs = 0;
		outBeams.push_back(outBeam);

		inBeam.id = newId;
		PushBeamToTree(inBeam, facetID, 0, Location::In);

#ifdef _CHECK_ENERGY_BALANCE
		ComputeFacetEnergy(facetID, outBeam);
#endif
	}

	TraceInternalBeams(outBeams);
}

void ScatteringConvex::ScatterLight(const std::vector<std::vector<int>> &/*tracks*/, std::vector<Beam> &)
{
<<<<<<< HEAD
}
=======
	m_incidentEnergy = 0;

	m_visibleFacets.nElems = 0;
	FindVisibleFacets(m_originalBeam, m_lightChecker, m_workFacets,
					  m_visibleFacets);
>>>>>>> origin/refactor

void ScatteringConvex::TraceInternalBeams(std::vector<Beam> &outBeams)
{
	while (m_treeSize != 0)
	{
<<<<<<< HEAD
		Beam beam = m_beamTree[--m_treeSize];

		if (IsTerminalAct(beam))
		{
			continue;
		}

		for (int id = 0; id < m_particle->nFacets; ++id)
		{
			if (id == beam.lastFacetId)
			{
				continue;
			}

			Beam inBeam;
			bool isIncident = SplitSecondaryBeams(beam, id, inBeam, outBeams);

			if (!isIncident)
			{
				continue;
			}

			inBeam.id = RecomputeTrackId(beam.id, id);
			inBeam.locations = beam.locations;
			PushBeamToTree(inBeam, id, beam.nActs+1, Location::In);
		}
=======
		Facet *facet = m_visibleFacets.elems[i];
		ComputeOpticalBeamParams(facet, m_originalBeam, *facet);
		PushBeamsToBuffer(facet, m_splitting.beams, externalBeams);
>>>>>>> origin/refactor
	}
}

bool ScatteringConvex::SplitSecondaryBeams(Beam &incidentBeam, int facetID,
										   Beam &inBeam, std::vector<Beam> &outBeams)
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

#ifdef _DEBUG // DEB
//			double p = m_splitting.ComputeOutgoingOpticalPath(outBeam);
//			outBeam.opticalPath += p;
//			outBeam.ops.push_back(p);
#else
#endif
			outBeam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
			outBeam.lastFacetId = facetID;
			outBeams.push_back(outBeam);
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
//#ifdef _DEBUG // DEB
//		double p = m_splitting.ComputeOutgoingOpticalPath(outBeam);
//		outBeam.opticalPath += p;
//		outBeam.ops.push_back(p);
//#else
//		outBeam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
//#endif
		outBeam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
		outBeam.lastFacetId = facetID;
		outBeams.push_back(outBeam);
	}

	return true;
}

double ScatteringConvex::MeasureOpticalPath(const Beam &beam,
											const Point3f sourcePoint,
											const std::vector<int> &track)
{
	double path = 0;
	Point3f dir = -beam.direction; // back direction
	Location loc = Location::Out;

	Point3f p1 = sourcePoint;
	Point3f p2;

	// back tracing
	for (int i = track.size()-1; i > 0; --i)
	{
		Point3f &exNormal = m_facets[track[i]].ex_normal;
		dir = m_splitting.ChangeBeamDirectionConvex(dir, exNormal, loc);

		Point3f &inNormal = m_facets[track[i-1]].in_normal;
		p2 = ProjectPointToPlane(p1, dir, inNormal);
		double len = Length(p2 - p1);

		path += len;
		p1 = p2;
		loc = Location::In;
	}

	lastPoint = p2;

	path *= real(m_splitting.GetRi());
	return path;
}

double ScatteringConvex::MeasureFullOpticalPath(const Beam &beam,
												const Point3f sourcePoint,
												const std::vector<int> &track)
{
	double path = MeasureOpticalPath(beam, sourcePoint, track);

	Point3f nFar1 = m_incidentDir;
	Point3f nFar2 = -beam.direction;
	double dd1 = m_splitting.FAR_ZONE_DISTANCE + DotProductD(lastPoint, nFar1);
	double dd2 = fabs(DotProductD(sourcePoint, nFar2) + m_splitting.FAR_ZONE_DISTANCE);
	path += dd1;
	path += dd2;

	return path;
=======
	if (m_hasTracks)
	{
		SetWorkFacetsFromTracks();
	}

	FindVisibleFacets(beam, m_lightChecker, m_workFacets, facets);
>>>>>>> origin/refactor
}
