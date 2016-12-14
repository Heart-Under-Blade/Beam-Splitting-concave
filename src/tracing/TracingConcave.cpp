#include "TracingConcave.h"

#include "macro.h"
#include <tgmath.h>
#include <assert.h>

#define MULTI_INDEX		10000000l			// index for poligon's clip operations
#define EPS_MULTI		(1.415*MULTI_INDEX*2)/10000	// погрешность, при которой точки операций Clipper'а можно считать совпадающими

#define SIMPLE_CLIP_RESULT 1

#define EPS_ORTO_FACET 0.0001 //TODO: подобрать норм значение

using namespace ClipperLib;

#ifdef _TRACK_ALLOW
std::ofstream trackMapFile("tracks.dat", std::ios::out);
#endif

TracingConcave::TracingConcave(Particle *particle, const Point3f &startBeamDir,
							   bool isOpticalPath, const Point3f &polarizationBasis,
							   int interReflectionNumber)
	: Tracing(particle, startBeamDir, isOpticalPath, polarizationBasis, interReflectionNumber)
{
	m_startBeam.polygon[0] = m_startBeam.direction * m_particle->halfHeight*2;
	++m_startBeam.size;
	m_startBeam.direction.d_param = DotProduct(m_startBeam.polygon[0], m_startBeam.direction);
	m_startBeam.isExternal = true;

	m_clipper.ZFillFunction(FindZCoord);
}

//bool TracingConcave::isOrderReversed(const Point3f oldNormal, const Path polygon)
//{
//	Point3f facet[MAX_VERTEX_NUM];

//	for (int i = 0; i < 3; ++i)
//	{
//		facet[i].cx = (float)polygon[i].X / MULTI_INDEX;
//		facet[i].cy = (float)polygon[i].Y / MULTI_INDEX;
//		facet[i].cz = (float)polygon[i].Z / MULTI_INDEX;
//	}

//	Point3f newNormal = NormalToFacet(facet);
//	double cosNO = DotProduct(newNormal, oldNormal);
//	return (cosNO < 0);
//}

void TracingConcave::InversePolygonOrder(Path &polygon)
{
	IntPoint tmp;
	int size = polygon.size()-1;

	for (int vi = 0; vi <= size/2; ++vi)
	{
		tmp = polygon[vi];
		polygon[vi] = polygon[size-vi];
		polygon[size-vi] = tmp;
	}
}

void TracingConcave::SplitBeamByParticle(std::vector<Beam> &outBeams,
										 double &lightSurfaceSquare)
{
	m_treeSize = 0;

	int facetIds[MAX_FACET_NUM];
	int facetIdNum = 0;

	SelectVisibleFacetsExternal(m_startBeam, facetIds, facetIdNum);
	SortFacets(facetIdNum, m_startBeam.direction, facetIds);

#ifdef _TRACK_OUTPUT
	trackMapFile << "0 lvl: ";
#endif

	// reflecting the beam to each visible facet
	for (int i = 0; i < facetIdNum; ++i)
	{
		int facetId = facetIds[i];

		Beam inBeam, outBeam;
		SetBeamsParamsExternal(facetId, inBeam, outBeam);

		if (i == 0) // this facet is not shadowed by others
		{
			SetBeamPolygonByFacet(facetId, inBeam);
			SetBeamPolygonByFacet(facetId, outBeam);
		}
		else // facet is probably shadowed by others
		{
			Paths cuttedFacet(1);

			CutShadowsFromFacet(facetId, facetIds, i, m_startBeam,
								cuttedFacet);

			if (!cuttedFacet.empty())
			{
				LOG_ASSERT(cuttedFacet.size() < 2);

				SetBeamByPath(inBeam, cuttedFacet.at(0));
				SetBeamByPath(outBeam, cuttedFacet.at(0));
			}
			else // facet is totaly shadowed by others
			{
				continue;
			}
		}

		if (m_isOpticalPath)
		{
			Point3f center = CenterOfPolygon(inBeam.polygon, inBeam.size);

			inBeam.D = DotProduct(-inBeam.direction, center);
			inBeam.opticalPath = FAR_ZONE_DISTANCE - DotProduct(-m_startBeam.direction, center);

			outBeam.D = DotProduct(-outBeam.direction, center);
			outBeam.opticalPath = inBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + outBeam.D);
		}

#ifdef _TRACK_OUTPUT
		trackMapFile << facetId << ", ";
#endif
		PushBeamToTree(inBeam, facetId, 0, false);
		PushBeamToTree(outBeam, facetId, 0, true);

//		lightSurfaceSquare += outBeam.Square() * cosIncident;
	}

#ifdef _TRACK_OUTPUT
	trackMapFile.flush();
#endif

	TraceInternalReflections(outBeams);
}

void TracingConcave::SelectVisibleFacets(const Beam &beam, int *facetIds, int &facetIdCount)
{
	FindVisibleFacetsInternal(beam, facetIds, facetIdCount);

	Point3f dir = beam.direction;
	dir.d_param = m_particle->normals[beam.facetId].d_param;
	SortFacets(facetIdCount, dir, facetIds);
}

void TracingConcave::CatchExternalBeam(const Beam &beam, std::vector<Beam> &outBeams)
{
	const Point3f &dir = beam.direction;
	Point3f &preNormal = m_particle->externalNormals[beam.facetId];

	int facetIds[MAX_FACET_NUM];
	int facetIdCount = 0;
	SelectVisibleFacets(beam, facetIds, facetIdCount);

	Paths originPolygon(1);
	SetPolygonByFacet(beam.polygon, beam.size, originPolygon);

	// cut facet projections out of beam one by one
	for (int i = 0; i < facetIdCount; ++i)
	{
		int id = facetIds[i];

		Paths clippedPolygon(1);
		CutBeamByFacet(originPolygon, id, dir, preNormal, clippedPolygon);

		if (clippedPolygon.empty()) /// beam incedents on facet totaly
		{
			originPolygon.clear();
			break;
		}
		else
		{
			originPolygon = clippedPolygon;
		}
	}

	Beam b = beam;

	for (const Path &p : originPolygon)
	{
		SetBeamByPath(b, p);
		outBeams.push_back(b);
	}
}

void TracingConcave::PushBeamToTree(Beam &beam, int facetId, int level,
									bool isExternal)
{
	assert(m_treeSize < MAX_BEAM_REFL_NUM);

#ifdef _TRACK_ALLOW
	std::vector<int> &tr = beam.track;
	int size = tr.size();

if (size == 0 || (size > 0 && facetId != tr.at(size-1)))
	tr.push_back(facetId);
#endif

#ifdef _TRACK_OUTPUT
	PrintTrack(beam, facetId);
#endif

	beam.facetId = facetId;
	beam.level = level;
	beam.isExternal = isExternal;
	m_beamTree[m_treeSize] = beam;
	++m_treeSize;
}

void TracingConcave::PrintTrack(const Beam &beam, int facetId)
{
#ifdef _TRACK_ALLOW
	trackMapFile << "(" << beam.track.at(0);

	for (int i = 1; i < beam.track.size(); ++i)
	{
		trackMapFile << ", ";
		trackMapFile << beam.track.at(i);
	}

	if (facetId != beam.track.back())
	{
		trackMapFile << ", " << facetId;
	}

	trackMapFile << "), ";
#endif
}

void TracingConcave::CutIncidentBeam(int facetId, Beam &beam, bool &isDivided)
{
	isDivided = false;

	const Point3f &previewNormal = (beam.isExternal)
			? m_particle->externalNormals[beam.facetId]
			: m_particle->normals[beam.facetId];

	Paths clippedFacet;
	Paths origin(1);
	SetPolygonByFacet(beam.polygon, beam.size, origin);
	CutBeamByFacet(origin, facetId, beam.direction, previewNormal, clippedFacet);

	if (clippedFacet.empty()) // beam is totaly swallowed by facet
	{
		beam.size = 0;
	}
	else if (clippedFacet.size() == SIMPLE_CLIP_RESULT)
	{
		SetBeamByPath(beam, clippedFacet.at(0));
	}
	else // beam had divided by facet
	{
		Beam b = beam;

		for (const Path &p : clippedFacet)
		{
			SetBeamByPath(b, p);
			PushBeamToTree(b, b.facetId, b.level, b.isExternal);
		}

		isDivided = true;
		beam.size = 0;
	}
}

void TracingConcave::TraceInternalReflections(std::vector<Beam> &outBeams)
{
	while (m_treeSize != 0)
	{
		Beam incidentBeam = m_beamTree[--m_treeSize];
		const bool isExternal = incidentBeam.isExternal;

		if (isTerminalBeam(incidentBeam))
		{
			if (isExternal) // отлов и отсечение отраженных пучков
			{
				CatchExternalBeam(incidentBeam, outBeams);
			}

			continue;
		}

#ifdef _TRACK_OUTPUT
		trackMapFile << "\n" << incidentBeam.level << " lvl: ";
		trackMapFile.flush();
#endif
		const Point3f &incidentDir = incidentBeam.direction;

		int facetIds[MAX_FACET_NUM];
		int facetIdCount = 0;
		SelectVisibleFacets(incidentBeam, facetIds, facetIdCount);

		for (int i = 0; i < facetIdCount; ++i)
		{
			int facetId = facetIds[i];

			Beam outBeam;
			bool isOk = Intersect(facetId, incidentBeam, outBeam);

			if (!isOk)
			{
				continue;
			}

			bool isDivided;
			CutIncidentBeam(facetId, incidentBeam, isDivided);

			Beam inBeam = outBeam;

#ifdef _TRACK_ALLOW
			outBeam.track = incidentBeam.track;
			inBeam.track = incidentBeam.track;
#endif
			const Point3f &normal = m_particle->externalNormals[facetId];
			double cosIN = DotProduct(incidentDir, normal);
			double cosI_sqr = cosIN * cosIN;

			double Nr;
			{
				double &re = m_particle->ri_coef_re;
				double &im = m_particle->ri_coef_im;
				Nr = (re + sqrt(re*re + im/cosI_sqr))/2.0;
			}

			if (cosIN >= EPS_COS_00) /// normal incidence
			{
				SetNormalIncidenceBeamParams(Nr, incidentBeam, inBeam, outBeam);
				PushBeamToTree(outBeam, facetId, incidentBeam.level+1, true);
			}
			else /// slopping incidence
			{
				if (!isExternal)
				{
					Point3f r0 = incidentDir/cosIN - normal;

					Point3f reflDir = r0 - normal;
					Normalize(reflDir);
					inBeam.direction = reflDir;

					Point3f scatteringNormal;
					CrossProduct(normal, incidentDir, scatteringNormal);
					Normalize(scatteringNormal);

					Beam incBeam = incidentBeam;
					incBeam.RotatePlane(scatteringNormal);

					double s = 1.0/(Nr*cosI_sqr) - Norm(r0);

					if (s > DBL_EPSILON)
					{
						SetTrivialIncidenceBeamParams(cosIN, Nr, normal, r0, s, incidentBeam,
													  inBeam, outBeam);
						PushBeamToTree(outBeam, facetId, incidentBeam.level+1, true);
					}
					else /// case of the complete internal reflection
					{
						SetCompleteReflectionBeamParams(cosIN, Nr, incidentBeam, inBeam);
					}
				}
				else // case of external beam incidents to facet
				{
					inBeam.JMatrix = incidentBeam.JMatrix;
					double cosI = DotProduct(-normal, incidentDir);

					SetSloppingBeamParamsExternal(incidentDir, cosI, facetId,
												  inBeam, outBeam);
					if (m_isOpticalPath)
					{
						CalcOpticalPathInternal(Nr, incidentBeam, inBeam, outBeam);
					}

					PushBeamToTree(outBeam, facetId, incidentBeam.level+1, true);
				}
			}

#ifdef _TRACK_ALLOW
			inBeam.track = incidentBeam.track;
#endif
			PushBeamToTree(inBeam, facetId, incidentBeam.level+1, false);

#ifdef _TRACK_OUTPUT
			trackMapFile << "[in], ";
#endif
			if (isDivided)
			{
				break;
			}
		}

		if (isExternal && incidentBeam.size != 0)
		{
			// посылаем обрезанный всеми гранями внешний пучок на сферу
			outBeams.push_back(incidentBeam);
		}
	}
}

// REF: объединить в одну со след.
void TracingConcave::SelectVisibleFacetsExternal(const Beam &beam, int *facetIndices,
												 int &indicesNumber)
{
	Point3f *normals = m_particle->normals;

	for (int i = 0; i < m_particle->facetNum; ++i)
	{
		double cosIncident = DotProduct(beam.direction, normals[i]);

		if (cosIncident > EPS_COS_90) /// beam incidents to this facet
		{
			facetIndices[indicesNumber] = i;
			++indicesNumber;
		}
	}
}

void TracingConcave::FindVisibleFacetsInternal(const Beam &beam, int *facetIndices,
											   int &facetIdCount)
{
	Point3f *normals = (beam.isExternal) ? m_particle->normals
										 : m_particle->externalNormals;

	const Point3f &invNormal = -normals[beam.facetId];

	Point3f cob = CenterOfPolygon(beam.polygon, beam.size);

	for (int i = 0; i < m_particle->facetNum; ++i)
	{
		double cosIncident = DotProduct(beam.direction, normals[i]);

		if (cosIncident >= EPS_COS_90) // beam incidents to this facet
		{
			Point3f cof = CenterOfPolygon(m_particle->facets[i].polygon, m_particle->facets[i].size);
			Point3f vectorToFacet = cof - cob;
			double cosFacets = DotProduct(invNormal, vectorToFacet);

			if (cosFacets >= EPS_ORTO_FACET) // facet is in front of begin of beam
			{
				facetIndices[facetIdCount] = i;
				++facetIdCount;
			}
		}
	}
}

void TracingConcave::SetPolygonByFacet(const Point3f *facet, int size, Paths &polygon) const
{
	for (int i = 0; i < size; ++i)
	{
		const Point3f &p = facet[i];

		polygon[0] << IntPoint((cInt)(p.cx * MULTI_INDEX),
							   (cInt)(p.cy * MULTI_INDEX),
							   (cInt)(p.cz * MULTI_INDEX));
	}
}

void TracingConcave::RemoveHole(Paths &result)
{
	Path &pol1 = result.front();
	Path &pol2 = result.back();

	bool or1 = ClipperLib::Orientation(pol1);
	bool or2 = ClipperLib::Orientation(pol2);

	if (or1 == or2)
	{
		return;
	}

	int first = 0, second = 0;
	cInt len;
	cInt min_len = LONG_MAX;
	IntPoint diff;

	for (int i = 0; i < pol1.size(); ++i)
	{
		for (int j = 0; j < pol2.size(); ++j)
		{
			diff.X = pol1.at(i).X - pol2.at(j).X;
			diff.Y = pol1.at(i).Y - pol2.at(j).Y;
			diff.Z = pol1.at(i).Z - pol2.at(j).Z;

			len = diff.X*diff.X + diff.Y*diff.Y + diff.Z*diff.Z;
			len = (cInt)std::round(sqrt((double)len));

			if (len < min_len)
			{
				first = i;
				second = j;
				min_len = len;
			}
		}
	}

	Paths newResult(1);

	for (int i = 0; i < pol1.size(); ++i)
	{
		newResult[0] << pol1.at(i);

		if (i == first)
		{
			int j = second;

			do
			{
				if (j == pol2.size())
				{
					j = 0;
				}
				else
				{
					newResult[0] << pol2.at(j);
					++j;
				}
			}
			while (j != second);

			newResult[0] << pol2.at(second);
			newResult[0] << pol1.at(first);
		}
	}

	result = newResult;
}

void TracingConcave::HandleResultPolygon(Axis axis, Paths &result)
{
	SwapCoords(Axis::aZ, axis, result); // обратно
	ClipperLib::CleanPolygons(result, EPS_MULTI);
	RemoveEmptyPolygons(result);

	LOG_ASSERT(result.size() < 3);

	if (result.size() == 2)
	{
		RemoveHole(result);
	}
}

void TracingConcave::CutShadowsFromFacet(int facetId, int *facetIds,
										 int previewFacetCount,
										 const Beam &incidentBeam,
										 Paths &resultPolygon)
{
	const Point3f *facet = m_particle->facets[facetId].polygon;
	const int size = m_particle->facets[facetId].size;

	const Point3f &originNormal = (incidentBeam.isExternal)
			? m_particle->externalNormals[facetId]
			  : m_particle->normals[facetId];

	Axis axis = GetSwapAxis(originNormal);

	SetPolygonByFacet(facet, size, resultPolygon); // set origin polygon
	SwapCoords(axis, Axis::aZ, resultPolygon);

	Paths clip(previewFacetCount);

	for (int i = 0; i < previewFacetCount; ++i)
	{
		// set clip polygon by projection
		int facetId = facetIds[i];
		const Point3f *fac = m_particle->facets[facetId].polygon;
		int size = m_particle->facets[facetId].size;
		ProjectFacetToFacet(fac, size, incidentBeam.direction, originNormal, clip[i]);

		{	// equate similar points /// REF: возможная замена ClearPolygon
//			for (IntPoint &p0 : resultPolygon[0])
//			{
//				for (IntPoint &p1 : clip[0])
//				{
//					if (abs(p1.X - p0.X) < EPS_MULTI
//						&& abs(p1.Y - p0.Y) < EPS_MULTI
//						&& abs(p1.Z - p0.Z) < EPS_MULTI)
//					{
//						p0 = p1;
//					}
//				}
//			}
		}
	}

	SwapCoords(axis, Axis::aZ, clip);
	ClipDifference(resultPolygon, clip, resultPolygon);

	RemoveEmptyPolygons(resultPolygon);

	if (!resultPolygon.empty())
	{
		HandleResultPolygon(axis, resultPolygon);
	}
}

void TracingConcave::ProjectPointToFacet(const Point3d &point, const Point3d &direction, const Point3d &facetNormal, Point3d &projection)
{
	double t = DotProductD(point, facetNormal);
	t = t + facetNormal.d;
	double dp = DotProductD(direction, facetNormal);
	t = t/dp;
	projection = point - (direction * t);
}

void TracingConcave::SetBeamByPath(Beam &beam, const Path &result)
{
	int vertexNum = result.size();
	beam.size = vertexNum;

	for (const IntPoint &p : result)
	{
		Point3f tmp((float)p.X / MULTI_INDEX,
					(float)p.Y / MULTI_INDEX,
					(float)p.Z / MULTI_INDEX);

		beam.polygon[--vertexNum] = tmp;
	}

	LOG_ASSERT(beam.size > 0 && beam.size <= MAX_VERTEX_NUM);
}

void TracingConcave::ProjectPointToFacet(const Point3f &point, const Point3f &direction,
										 const Point3f &facetNormal, Point3f &projection)
{
	double t = DotProduct(point, facetNormal);
	t = t + facetNormal.d_param;
	double dp = DotProduct(direction, facetNormal);
	t = t/dp;
	projection = point - (direction * t);
}

/// OPT: поменять все int и пр. параметры функций на ссылочные

void TracingConcave::ProjectFacetToFacet(const Point3f *a_facet, int a_size,
										 const Point3f &a_dir,
										 const Point3f &b_normal,
										 Path &projection)
{
	for (int i = 0; i < a_size; ++i)
	{
		Point3d p;
		ProjectPointToFacet(Point3d(a_facet[i].cx, a_facet[i].cy, a_facet[i].cz),
							Point3d(a_dir.cx, a_dir.cy, a_dir.cz),
							Point3d(b_normal.cx, b_normal.cy, b_normal.cz, b_normal.d_param), p);

		projection << IntPoint((cInt)(p.x * MULTI_INDEX),
							   (cInt)(p.y * MULTI_INDEX),
							   (cInt)(p.z * MULTI_INDEX));
	}
}

double TracingConcave::CalcMinDistanceToFacet(int facetId, const Point3f &beamDir)
{
	double dist = FLT_MAX;
	Point3f *facet = m_particle->facets[facetId].polygon;

	for (int i = 0; i < m_particle->facets[facetId].size; ++i)
	{
		// measure dist
		Point3f point;
		ProjectPointToFacet(facet[i], -beamDir, beamDir, point);
		double newDist = sqrt(Norm(point - facet[i]));

		if (newDist < dist) // choose minimum with previews
		{
			dist = newDist;
		}
	}

	return dist;
}

void TracingConcave::SortFacets(int number, const Point3f &beamDir, int *facetIds)
{
	float distances[MAX_VERTEX_NUM];

	for (int i = 0; i < number; ++i)
	{
		distances[i] = CalcMinDistanceToFacet(facetIds[i], beamDir);
	}

	int left = 0;
	int rigth = number - 1;

	int stack[MAX_VERTEX_NUM*2];
	int size = 0;

	stack[size++] = left;
	stack[size++] = rigth;

	while (true)
	{
		float base = distances[(left + rigth)/2];

		int i = left;
		int j = rigth;

		while (i <= j)
		{
			while (distances[i] < base)
			{
				++i;
			}

			while (distances[j] > base)
			{
				--j;
			}

			if (i <= j)	// exchange elems
			{
				float temp_d = distances[i];
				distances[i] = distances[j];
				distances[j] = temp_d;

				int temp_v = facetIds[i];
				facetIds[i] = facetIds[j];
				facetIds[j] = temp_v;

				++i;
				--j;
			}
		}

		if (i < rigth)
		{
			stack[size++] = i;
			stack[size++] = rigth;
		}

		if (left < j)
		{
			stack[size++] = left;
			stack[size++] = j;
		}

		if (size == 0)
		{
			break;
		}

		rigth = stack[--size];
		left = stack[--size];
	}
}

void TracingConcave::SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
										 std::vector<Beam> &/*outBeams*/)
{
	for (unsigned int i = 0; i < tracks.size(); ++i)
	{
		/// TODO: realize
	}
}

void TracingConcave::SwapCoords(Axis oldAxis, Axis newAxis, Paths &origin) const
{
	if (oldAxis == newAxis)
	{
		return;
	}

	cInt *oldP, *newP;

	for (Path &path : origin)
	{
		for (IntPoint &p : path)
		{
			// REF: do smth
			switch (oldAxis) {
			case Axis::aX: oldP = &p.X;
				break;
			case Axis::aY: oldP = &p.Y;
				break;
			case Axis::aZ: oldP = &p.Z;
				break;
			}

			switch (newAxis) {
			case Axis::aX: newP = &p.X;
				break;
			case Axis::aY: newP = &p.Y;
				break;
			case Axis::aZ: newP = &p.Z;
				break;
			}

			cInt oldP1 = *oldP;
			cInt newP1 = *newP;
			*oldP = newP1;
			*newP = oldP1;
		}
	}
}

Axis TracingConcave::GetSwapAxis(const Point3f &normal)
{
	if (fabs(normal.cx) > 0.5)
	{
		return Axis::aX;
	}
	else if (fabs(normal.cy) > 0.5)
	{
		return Axis::aY;
	}
	else // no need to swap
	{
		return Axis::aZ;
	}
}

void TracingConcave::ClipDifference(const Paths &subject, const Paths &clip,
									Paths &difference)
{
	m_clipper.AddPaths(subject, ptSubject, true);
	m_clipper.AddPaths(clip, ptClip, true);
	m_clipper.Execute(ctDifference, difference);
	m_clipper.Clear();
}

void TracingConcave::CutBeamByFacet(Paths &beamPolygon, int facetId,
									const Point3f &direction,
									const Point3f &shapeNormal,
									Paths &result)
{
	Axis axis = GetSwapAxis(shapeNormal);
	SwapCoords(axis, Axis::aZ, beamPolygon);

	Paths clip(1);
	{
		const Point3f *facet = m_particle->facets[facetId].polygon;
		int size = m_particle->facets[facetId].size;
		ProjectFacetToFacet(facet, size, direction, shapeNormal, clip[0]); // проецируем грань на начальный пучок
		SwapCoords(axis, Axis::aZ, clip);
	}

	ClipDifference(beamPolygon, clip, result);

	if (!result.empty())
	{
		HandleResultPolygon(axis, result);
	}
}

void TracingConcave::RemoveEmptyPolygons(Paths &result)
{
	Paths buff = result;
	result.clear();

	for (int i = 0; i < buff.size(); ++i)
	{
		if (buff.at(i).size() >= SIMPLE_CLIP_RESULT)
		{
			result.push_back(buff.at(i));
		}
	}
}

void TracingConcave::FindDividePoint(const std::vector<Point3f> &polygon,
									 int i0, int i1, const Point3f &normal,
									 Point3f &x, int &nextPointIndex) const
{
	int size = polygon.size();
	__m128 res;

	__m128 _a0 = _mm_setr_ps(polygon[i0].cx, polygon[i0].cy, polygon[i0].cz, 0.0);
	__m128 _a1 = _mm_setr_ps(polygon[i1].cx, polygon[i1].cy, polygon[i1].cz, 0.0);
	__m128 _n = _mm_load_ps(normal.point);

	bool isOk = false;
	int j = size-1;

	for (int i = 0; i < size; ++i)
	{
		__m128 _b0 = _mm_setr_ps(polygon[j].cx, polygon[j].cy, polygon[j].cz, 0.0);
		__m128 _b1 = _mm_setr_ps(polygon[i].cx, polygon[i].cy, polygon[i].cz, 0.0);

		if (is_inside_i(_b1, _a0, _a1, _n)
				&& !is_inside_i(_b0, _a0, _a1, _n))
		{
			res = computeIntersection_i(_a0, _a1, _b0, _b1, _n, isOk); // OPT: написать опт. вариант (с уже готовыми векторами v1, v2)

			if (isOk)
			{
				x.point[0] = res[0];
				x.point[1] = res[1];
				x.point[2] = res[2];
				nextPointIndex = i;
				return;
			}
		}

		j = i;
	}

	LOG_ASSERT(false && "Divide point is not found");
}

void TracingConcave::FillSubpolygon(int begin, int end,
									const std::vector<Point3f> &polygon,
									std::vector<Point3f> &subpolygon) const
{
	for (int j = begin; j != end; ++j)
	{
		if (j == polygon.size())
		{
			j = -1;
			continue;
		}

		subpolygon.push_back(polygon[j]);
	}
}

void TracingConcave::DividePolygon(const std::vector<Point3f> &polygon,
								   const Point3f &normal, Polygons &polygons) const
{
	int size = polygon.size();
	int baseI = size-1;

	for (int i = 1; i < size; ++i)
	{
		Point3f v1 = polygon[i-1] - polygon[baseI];
		Point3f v2 = polygon[i] - polygon[baseI];
		Point3f res;
		CrossProduct(v1, v2, res);
		double dir = DotProduct(res, normal);

		if (dir < -EPS_INTERSECTION) // cavity is here
		{
			Point3f x;
			int next;
			FindDividePoint(polygon, baseI, i-1, normal,
							x, next);

			auto Divide = [&] (int begin, int end)
			{
				std::vector<Point3f> subpolygon;
				subpolygon.push_back(polygon[end]);
				subpolygon.push_back(x);
				FillSubpolygon(begin, end, polygon, subpolygon);
				DividePolygon(subpolygon, normal, polygons);
			};

			Divide(next, i-1);
			Divide(i-1, next);

			return;
		}

		baseI = i-1;
	}

	polygons.push_back(polygon);
}

void TracingConcave::DivideConcavePolygon(const Point3f *polygon, int size,
										  const Point3f &normal,
										  Polygons &polygons) const
{
	std::vector<Point3f> vec_polygon;

	for (int i = 0; i < size; ++i)
	{
		vec_polygon.push_back(polygon[i]);
	}

	DividePolygon(vec_polygon, normal, polygons); // recursive
}

double TracingConcave::SquareOfPolygon(const std::vector<Point3f> &polygon) const
{
	double square = 0;
	const Point3f &basePoint = polygon[0];
	Point3f v1 = polygon[1] - basePoint;

	for (int i = 2; i < polygon.size(); ++i)
	{
		Point3f v2 = polygon[i] - basePoint;
		Point3f res;
		CrossProduct(v1, v2, res);
		square += sqrt(Norm(res));
		v1 = v2;
	}

	return square/2;
}

void FindZCoord(IntPoint &a1, IntPoint &a2, IntPoint &, IntPoint &, IntPoint &point)
{
	IntPoint *top, *bot;

	if (a2.Z > a1.Z)
	{
		top = &a2;
		bot = &a1;
	}
	else
	{
		top = &a1;
		bot = &a2;
	}

	double x, y;

	x = top->X - bot->X;
	y = top->Y - bot->Y;
	double normVec = x*x + y*y;

	x = point.X - bot->X;
	y = point.Y - bot->Y;
	double normBotVec = x*x + y*y;

	point.Z = (top->Z - bot->Z) * sqrt(normBotVec) / sqrt(normVec) + bot->Z;
}

double TracingConcave::AreaByClipper(const Beam &beam, const Point3f &normal) const
{
	double area = 0;

	Paths polygon(1);
	SetPolygonByFacet(beam.polygon, beam.size, polygon);

	Point3f normal1 = normal;

	Point3f n_normal = normal1;
	Normalize(n_normal); // REF: перегрузить функцию, чтобы возвращала значение

	if (fabs(n_normal.cx) > 0.5)
	{
		SwapCoords(Axis::aX, Axis::aZ, polygon);
		float tmp = normal1.cx;
		normal1.cx = normal1.cz;
		normal1.cz = tmp;
	}
	else if (fabs(n_normal.cy) > 0.5)
	{
		SwapCoords(Axis::aY, Axis::aZ, polygon);
		float tmp = normal1.cy;
		normal1.cy = normal1.cz;
		normal1.cz = tmp;
	}

	area = ClipperLib::Area(polygon[0]);
	area /= (cInt)MULTI_INDEX * MULTI_INDEX;

	Point3f z(0, 0, 1);
	double cosNormZ = DotProduct(normal1, z);

	if (cosNormZ > 0)
	{
		cosNormZ = DotProduct(normal1, -z);
	}

	area = (area*sqrt(Norm(normal1)))/cosNormZ;
//	Polygons polygons;
//	DivideConcavePolygon(beam.shape, beam.shapeSize, normal,
//						 polygons);

//	for (unsigned int i = 0; i < polygons.size(); ++i)
//	{
//		square += SquareOfPolygon(polygons[i]);
//	}

	if (area < 0)
	{
		area = -area;
	}

	return area;
}

double TracingConcave::BeamCrossSection(const Beam &beam) const
{
	const double eps = 1e7*DBL_EPSILON;

//	Point3f normal;
//	Point3f p1 = beam.polygon[1] - beam.polygon[0];
//	Point3f p2 = beam.polygon[2] - beam.polygon[0];
//	CrossProduct(p1, p2, normal);

	// normal of last facet of beam
	Point3f normal = m_particle->externalNormals[beam.facetId];

	double cosND = DotProduct(normal, beam.direction);
//	Normalize(normal);

//	if (fabs(normal.cz) < 0.5)
//	{
//		int minI = 0;

//		for (int i = 1; i < beam.shapeSize; ++i)
//		{
//			if (beam.shape[i].cz < beam.shape[minI].cz)
//			{
//				minI = i;
//			}
//		}
//	}
//	else if ...

//	if (dir < -FLT_EPSILON) /// нормаль посчитана неверно, разворачиваем её
//	{
//		normal.cx = -normal.cx;
//		normal.cy = -normal.cy;
//		normal.cz = -normal.cz;
//	}

	double e = fabs(cosND);

	if (e < eps)
	{
		return 0;
	}

	double square = AreaByClipper(beam, normal);
	double n = Length(normal);
	return (e*square) / n;
}
