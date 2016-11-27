#include "TracingConcave.h"

#include <assert.h>
#include <tgmath.h>
#include <iostream> // DEB
#include <fstream> // DEB

#define MULTI_INDEX		10000000l			// index for poligon's clip operations
#define EPS_MULTI		(1.415*MULTI_INDEX*2)/10000	// погрешность, при которой точки операций Clipper'а можно считать совпадающими

#define MAX_POLYGON_RESULT 1

using namespace ClipperLib;

static int outBeamNum = 0; // REF: del

std::ofstream trackMapFile("WW.dat", std::ios::out);

TracingConcave::TracingConcave(Particle *particle, const Point3f &startBeamDir,
							   bool isOpticalPath, const Point3f &polarizationBasis,
							   int interReflectionNumber)
	: Tracing(particle, startBeamDir, isOpticalPath, polarizationBasis, interReflectionNumber)
{
	m_startBeam.polygon[0] = -m_startBeam.direction * m_particle->halfHeight*2;
	++m_startBeam.size;
	m_startBeam.direction.d_param = DotProduct(m_startBeam.polygon[0], m_startBeam.direction);
	m_startBeam.isExternal = true;

	m_clipper.ZFillFunction(MeasureZ);
}

bool TracingConcave::isOrderReversed(const Point3f oldNormal, const Path polygon)
{
	Point3f facet[MAX_VERTEX_NUM];

	for (int i = 0; i < 3; ++i)
	{
		facet[i].cx = (float)polygon[i].X / MULTI_INDEX;
		facet[i].cy = (float)polygon[i].Y / MULTI_INDEX;
		facet[i].cz = (float)polygon[i].Z / MULTI_INDEX;
	}

	Point3f newNormal = NormalToFacet(facet);
	double cosNO = DotProduct(newNormal, oldNormal);
	return (cosNO < 0);
}

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

void TracingConcave::PrepareVisibleFacets(const Beam &beam, int *facetIds,
										  int &idNumber)
{
//	FindVisibleFacets(beam, facetIds, idNumber);
//	SortFacets(idNumber, beam.direction, facetIds);
}

void TracingConcave::SplitBeamByParticle(std::vector<Beam> &outBeams,
										 double &lightSurfaceSquare)
{
	m_treeSize = 0;

	int facetIds[MAX_FACET_NUM];
	int idNum = 0;
//	PrepareVisibleFacets(m_startBeam, facetIds, idNum);
	SelectVisibleFacetsExternal(m_startBeam, facetIds, idNum);
	SortFacets(idNum, m_startBeam.direction, facetIds);

#ifdef _OUTPUT_TRACK
	trackMapFile << "0 lvl: "; // DEB
#endif

	// reflecting the beam to each visible facet
	for (int i = 0; i < idNum; ++i)
	{
		int facetId = facetIds[i];
		const Point3f &extNormal = m_particle->externalNormals[facetId];
		double cosBN = DotProduct(m_startBeam.direction, extNormal);

		Beam inBeam, outBeam;
		SetBeamsParamsExternal(facetId, cosBN, inBeam, outBeam);

		if (i == 0) // this facet is not shadowed by others
		{
			SetBeamByFacet(facetId, inBeam);
			SetBeamByFacet(facetId, outBeam);
		}
		else // facet is probably shadowed by others
		{
			Paths cuttedFacet(1);

			const Point3f *facet = m_particle->facets[facetId];
			int size = m_particle->vertexNums[facetId];
			m_startBeam.facetId = facetId;
			CutShadowsFromFacet(facet, size, facetIds, i, m_startBeam,
								cuttedFacet);

			if (cuttedFacet.empty() //||
					/*WTF*//*(!cuttedFacet.empty() && cuttedFacet.at(0).empty())*/) /// facet is totaly shadowed by others
			{
				continue;
			}

			if (cuttedFacet.size() >= 2)
			assert(cuttedFacet.size() < 2);

			if (isOrderReversed(extNormal, cuttedFacet.at(0)))
			{
				InversePolygonOrder(cuttedFacet.at(0));
			}

			SetBeamShapeByPolygon(inBeam, cuttedFacet.at(0));
			SetBeamShapeByPolygon(outBeam, cuttedFacet.at(0));
		}

		if (m_isOpticalPath)
		{
			Point3f center = CenterOfPolygon(inBeam.polygon, inBeam.size);

			inBeam.D = DotProduct(-inBeam.direction, center);
			inBeam.opticalPath = FAR_ZONE_DISTANCE - DotProduct(m_startBeam.direction, center);

			outBeam.D = DotProduct(-outBeam.direction, center);
			outBeam.opticalPath = inBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + outBeam.D);
		}

#ifdef _OUTPUT_TRACK
		trackMapFile << facetId << ", ";
#endif
		PushBeamToTree(inBeam, facetId, 0, false);
		PushBeamToTree(outBeam, facetId, 0, true);
//		lightSurfaceSquare += outBeam.Square() * cosIncident;
	}

#ifdef _OUTPUT_TRACK
	trackMapFile.flush();
#endif

	TraceInternalReflections(outBeams);

	//DEB
//	int co = 0;
//	double s = 0;
//	double d = 0;
//	for (const Beam &b : outBeams)
//	{
//		switch (b.track[0]) {
//		case 0: case 7:
//			if (b.track.size() > 1 /*&& b.track[1] == 5*/)
//				s += Square(b);
//			break;
////		case 2:
//		default:
//			if (b.track.size() > 2)
//				d += Square(b);
//		}

//		if (b.track.size() == 3)
//		{
//			for (int a : b.track)
//				std::cout << a;
//		std::cout << std::endl;
//			co++;
//			d += Square(b);
//			std::cout << d << std::endl;
//		}
//	}

//	int fff = 0;
}

// REF: del
void TracingConcave::CutReflectedBeam(const Beam &beam, Beam &incidentBeam)
{
//	const Point3f &incidentDir = beam.direction;
//	Point3f &previewNormal = m_particle->externalNormals[beam.facetId];

//	int facetIds[MAX_FACET_NUM]; // OPT: заменить макс размер на 18
//	int idCount = 0;
//	SelectVisibleFacetsInternal(beam, facetIds, idCount);

//	incidentDir.D_PARAM = m_particle->normals[beam.facetId].D_PARAM;
//	SortFacets(idCount, incidentDir, facetIds);

//	Paths origin(1);
//	SetPolygonByFacet(beam.polygon, beam.size, origin);

//	for (int i = 0; i < idCount; ++i)
//	{
//		int facetId = facetIds[i];

//		Paths cuttedFacet(1);
//		CutBeamShapeByFacet(beam, facetId, incidentDir, previewNormal, cuttedFacet);

//		if (cuttedFacet.empty()) /// beam is totaly swallowed by facet
//		{
//			break;
//		}
//		else/*if (cuttedFacet.size() == MAX_POLYGON_RESULT
//				 && cuttedFacet.at(0).size() >= MIN_VERTEX_NUM)*/
//		{
//			origin = cuttedFacet;
//		}
//	}

//	for (const Path &p : origin)
//	{
//		Beam b;
//		SetBeamShapeByPolygon(b, p);
//		outBeams.push_back(b);
//	}
}

void TracingConcave::CatchExternalBeam(const Beam &beam, std::vector<Beam> &outBeams)
{
	Point3f incidentDir = beam.direction; // OPT: ссылку
	Point3f &previewNormal = m_particle->externalNormals[beam.facetId];

	int facetIds[MAX_FACET_NUM]; // OPT: заменить макс размер на 18
	int idCount = 0;
	SelectVisibleFacetsInternal(beam, facetIds, idCount);

	incidentDir.d_param = m_particle->normals[beam.facetId].d_param;
	SortFacets(idCount, incidentDir, facetIds);

	Paths origin(1);
	SetPolygonByFacet(beam.polygon, beam.size, origin);

	if (beam.track.size() == 5
			&& beam.track.at(4) == 12
			&& beam.track.at(3) == 11
			&& beam.track.at(2) == 8
			&& beam.track.at(1) == 11
			&& beam.track.at(0) == 8)
		int fff = 0;
	// cut facet projections out of beam one by one
	for (int i = 0; i < idCount; ++i)
	{
		int facetId = facetIds[i];

		Paths clippedFacet(1);
		CutBeamShapeByFacet(origin, facetId, incidentDir, previewNormal, clippedFacet);

		if (clippedFacet.empty()) /// beam is totaly swallowed by facet
		{
			break;
		}
		else/*if (cuttedFacet.size() == MAX_POLYGON_RESULT
				 && cuttedFacet.at(0).size() >= MIN_VERTEX_NUM)*/
		{
			origin = clippedFacet;
		}
	}

//	if (origin.size() >= 2)
	assert(origin.size() < 2);

	for (const Path &p : origin)
	{
		Beam b = beam;
		SetBeamShapeByPolygon(b, p);
		outBeams.push_back(b);
	}

//	if (incidentBeam.size != 0) // посылаем обрезанный всеми гранями внешний пучок на сферу
//	{
//		outBeams.push_back(incidentBeam);
//		++outBeamNum;//DEB
//	}
}

void TracingConcave::PushBeamToTree(Beam &beam, int facetId, int level,
									bool isExternal)
{
#ifdef _WRITE_TRACK
	std::vector<int> &tr = beam.track;
	int size = tr.size();

if (size == 0 || (size > 0 && facetId != tr.at(size-1)))
	tr.push_back(facetId);
#endif

#ifdef _OUTPUT_TRACK
	printTrack(beam, facetId);
#endif

	beam.facetId = facetId;
	beam.level = level;
	beam.isExternal = isExternal;
	m_tree[m_treeSize] = beam;
	++m_treeSize;
}

void TracingConcave::printTrack(const Beam &beam, int facetId)
{
#ifdef _WRITE_TRACK
	trackMapFile << "(" << beam.track.at(0);

	for (int i = 1; i < beam.track.size(); ++i)
	{
		trackMapFile << ", ";
		trackMapFile << beam.track.at(i);
	}

	if (facetId != beam.track.back())
	{
		trackMapFile << ", " << facetId << "), ";
	}
#endif
}

void TracingConcave::PushOutputBeamToTree(Beam &outBeam, std::vector<Beam> &buff,
										  int facetId, bool isDivided,
										  const Beam &incidentBeam)
{
	if (isDivided)
	{
		for (Beam &b : buff)
		{
			m_tree[m_treeSize++] = b;
		}

		buff.clear();
	}
	else
	{
#ifdef _WRITE_TRACK
		outBeam.track = incidentBeam.track;
#endif
		PushBeamToTree(outBeam, facetId, incidentBeam.level+1, true);
	}
//	PushBeamToTree(outBeam, facetId, incidentBeam.dept+1, true);
}

void TracingConcave::TraceInternalReflections(std::vector<Beam> &outBeams)
{
	std::vector<Beam> buff; // для хранения разделённых пучков

	while (m_treeSize != 0)
	{
		assert(m_treeSize < MAX_BEAM_REFL_NUM);

		Beam incidentBeam = m_tree[--m_treeSize]; // OPT: попробовать поменять на ссылку

		if (isEnough(incidentBeam))
		{
			if (incidentBeam.isExternal) // отлов и отсечение отраженных пучков (для зеркалки)
			{
				CatchExternalBeam(incidentBeam, outBeams);
			}

			continue;
		}

//if (incidentBeam.track.size() >= 2
//		&& incidentBeam.track.at(3) == 17
//		&& incidentBeam.track.at(2) == 9
//		&& incidentBeam.track.at(1) == 1
//		&& incidentBeam.track.at(0) == 1) //DEB
//	int fff = 0;

		Point3f &incidentDir = incidentBeam.direction;
		const bool isExternal = incidentBeam.isExternal;

		int facetIds[MAX_FACET_NUM]; // OPT: заменить макс размер на 18
		int idCount = 0;
		incidentDir.d_param = m_particle->normals[incidentBeam.facetId].d_param;

		// OPT: выполнять по условию только если есть видимые вершины
//		PrepareVisibleFacets(beam, facetIds, idCount);
		SelectVisibleFacetsInternal(incidentBeam, facetIds, idCount);
		SortFacets(idCount, incidentDir, facetIds);

#ifdef _OUTPUT_TRACK
		trackMapFile << "\n" << incidentBeam.level << " lvl: ";//DEB
		trackMapFile.flush();
#endif
		for (int i = 0; i < idCount; ++i)
		{
			int facetId = facetIds[i];

			if (facetId == incidentBeam.facetId) // same facet
			{
				continue; // OPT: проверить, проходит ли хоть раз эта проверка
			}

			Beam inBeam, outBeam;

			const Point3f &normal = m_particle->externalNormals[facetId];
			double cosIncident = DotProduct(normal, incidentDir);

			bool isOk = Intersect(facetId, incidentBeam, outBeam);

			if (!isOk)
			{
				continue;
			}

			bool isDivided = false;

			if (i != 0)
			{
				Paths clippedBeam(1);
				Beam bi = incidentBeam;
				bi.isExternal = true;
				bi.facetId = facetId;

				// обрезаем пучок попадающий на грань facetId
				CutShadowsFromFacet(outBeam.polygon, outBeam.size, facetIds,
									i, bi, clippedBeam);

				if (clippedBeam.empty()) /// facet is totaly shadowed by others (beam do not incedent on it)
				{
					continue;
				}
				else
				{
					if (clippedBeam.size() == MAX_POLYGON_RESULT)
					{
						if (isOrderReversed(NormalToFacet(outBeam.polygon), clippedBeam.at(0)))
						{// OPT: попробовать удалить (возможно лишнее)
							InversePolygonOrder(clippedBeam.at(0));
						}

						SetBeamShapeByPolygon(outBeam, clippedBeam.at(0));
					}
					else // beam had divided in several parts by facet
					{
						isDivided = true;

						for (const Path &p : clippedBeam)
						{
							Beam b = incidentBeam;
							SetBeamShapeByPolygon(b, p);
							buff.push_back(b);
						}

						incidentBeam.size = 0;
					}
//					if (isOrderReversed(NormalToFacet(outBeam.polygon), clippedBeam.at(0)))
//					{// OPT: попробовать удалить (возможно лишнее)
//						InversePolygonOrder(clippedBeam.at(0));
//					}

//					SetBeamShapeByPolygon(outBeam, clippedBeam.at(0));
				}
			}

			bool isIncidentDivided = false;

			if (isExternal)
			{
				const Point3f &previewNormal = m_particle->externalNormals[incidentBeam.facetId];
				Paths clippedFacet;
				CutBeamShapeByFacet(facetId, incidentBeam, previewNormal, clippedFacet);

				if (clippedFacet.empty()) /// beam is totaly swallowed by facet
				{
					incidentBeam.size = 0;
				}
				else /*if (cuttedFacet.size() == MAX_POLYGON_RESULT
						 && cuttedFacet.at(0).size() >= MIN_VERTEX_NUM)*/
				{
					if (clippedFacet.size() == MAX_POLYGON_RESULT)
					{
						if (isOrderReversed(NormalToFacet(incidentBeam.polygon), clippedFacet.at(0)))
						{// OPT: попробовать удалить (возможно лишнее)
							InversePolygonOrder(clippedFacet.at(0));
						}

						SetBeamShapeByPolygon(incidentBeam, clippedFacet.at(0));
					}
					else // beam had divided by facet
					{
						for (const Path &p : clippedFacet)
						{
							Beam b = incidentBeam;
							SetBeamShapeByPolygon(b, p);
							PushBeamToTree(b, b.facetId, b.level, b.isExternal);
//							outBeams.push_back(b);
						}

						isIncidentDivided = true;
						incidentBeam.size = 0;
					}
				}
			}

			inBeam = outBeam;

			double cosI_sqr = cosIncident * cosIncident;

			double Nr;
			{
				double &re = m_particle->refrI_coef_re;
				double &im = m_particle->refrI_coef_im;
				Nr = (re + sqrt(re*re + im/cosI_sqr))/2.0;
			}

			const complex &refrIndex = m_particle->refractionIndex;

			if (cosIncident >= EPS_COS_00) /// normal incidence
			{
				complex temp;

				temp = (2.0*refrIndex)/(1.0 + refrIndex); /// OPT: выделить эту константу
				SetBeam(outBeam, incidentBeam, incidentDir, incidentBeam.e,
						temp, temp);

				temp = (1.0 - refrIndex)/(1.0 + refrIndex);
				SetBeam(inBeam, incidentBeam, -incidentDir, incidentBeam.e,
						temp, -temp);

				if (m_isOpticalPath)
				{
					CalcOpticalPathInternal(Nr, incidentBeam, inBeam, outBeam);
				}

				PushOutputBeamToTree(outBeam, buff, facetId, isDivided, incidentBeam);
			}
			else /// slopping incidence
			{
				if (!isExternal)
				{
					Point3f r0 = incidentDir/cosIncident - normal;

					Point3f reflDir = r0 - normal;
					Normalize(reflDir);
					inBeam.direction = reflDir;

					Point3f scatteringNormal;
					CrossProduct(normal, incidentDir, scatteringNormal);
					Normalize(scatteringNormal);

					Beam incBeam = incidentBeam;
					incBeam.RotatePlane(scatteringNormal);

					inBeam.e = scatteringNormal;

					double s = 1.0/(Nr*cosI_sqr) - Norm(r0);
					complex tmp0 = refrIndex*cosIncident;

					if (s > DBL_EPSILON)
					{
						Point3f refrDir = r0/sqrt(s) + normal;
						Normalize(refrDir);

						double cosRefr = DotProduct(normal, refrDir);

						complex tmp1 = refrIndex*cosRefr;
						complex tmp = 2.0*tmp0;
						complex Tv0 = tmp1 + cosIncident;
						complex Th0 = tmp0 + cosRefr;

						SetBeam(outBeam, incBeam, refrDir, scatteringNormal,
								tmp/Tv0, tmp/Th0);

						complex Tv = (cosIncident - tmp1)/Tv0;
						complex Th = (tmp0 - cosRefr)/Th0;
						SetBeam(inBeam, incBeam, reflDir, scatteringNormal, Tv, Th);

						if (m_isOpticalPath)
						{
							CalcOpticalPathInternal(Nr, incBeam, inBeam, outBeam);
						}

						PushOutputBeamToTree(outBeam, buff, facetId, isDivided, incidentBeam);
					}
					else /// case of the complete internal reflection
					{
						const double bf = Nr*(1.0 - cosI_sqr) - 1.0;

						double im = (bf > 0) ? sqrt(bf) : 0;

						const complex sq(0, im);
						complex tmp = refrIndex*sq;
						complex Rv = (cosIncident - tmp)/(tmp + cosIncident);
						complex Rh = (tmp0 - sq)/(tmp0 + sq);

						SetBeam(inBeam, incidentBeam, reflDir, scatteringNormal, Rv, Rh);

						if (m_isOpticalPath)
						{
							Point3f center = CenterOfPolygon(outBeam.polygon, outBeam.size);
							inBeam.D = DotProduct(-center, inBeam.direction);

							double temp = DotProduct(incidentDir, center);

							inBeam.opticalPath = incidentBeam.opticalPath
									+ sqrt(Nr)*fabs(temp + incidentBeam.D);
						}
					}
				}
				else /// case of external beam incidents to facet
				{
//					/// rotate polarization plane // TODO: доделать
//					{
//						Point3f newBasis;
//						CrossProduct(normal, -incidentDir, newBasis);
//						Normalize(newBasis);
//						inBeam.e = ???;
//						inBeam.direction = incidentDir;
//						inBeam.RotatePlane(newBasis);
//					}

					Point3f refrDir, reflDir;

//					double cosI = cosIncident;
//					Point3f incDir = incidentDir;

					double cosI = DotProduct(normal, -incidentDir);// NOTE: используется косинус между двумя сонаправленными векторами (в отличие от остальных случаев)
					Point3f incDir = -incidentDir;

					SplitBeamDirection(incDir, cosI, normal,
									   reflDir, refrDir);

					double cosRelr = DotProduct(normal, reflDir);

					complex Tv00 = refrIndex*cosIncident;
					complex Th00 = refrIndex*cosRelr;

					complex Tv0 = Tv00 + cosRelr;
					complex Th0 = Th00 + cosIncident;

					complex Tv = (Tv00 - cosRelr)/Tv0; // OPT
					complex Th = (cosIncident - Th00)/Th0;
					SetBeam(outBeam, inBeam, refrDir, inBeam.e, Tv, Th);

					double cosInc2 = (2.0*cosIncident);
					SetBeam(inBeam, inBeam, reflDir, inBeam.e,
							cosInc2/Tv0, cosInc2/Th0);

					if (m_isOpticalPath)
					{
						CalcOpticalPathInternal(Nr, incidentBeam, inBeam, outBeam);
					}

					PushOutputBeamToTree(outBeam, buff, facetId, isDivided, incidentBeam);
				}
			}

			if (isDivided)
			{
				continue;
			}

#ifdef _WRITE_TRACK
			inBeam.track = incidentBeam.track;
#endif
			PushBeamToTree(inBeam, facetId, incidentBeam.level+1, false);

#ifdef _OUTPUT_TRACK
			trackMapFile << "[in], ";
#endif
			if (isIncidentDivided)
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
	Point3f *normals = m_particle->externalNormals;

	for (int i = 0; i < m_particle->facetNum; ++i)
	{
		double cosIncident = DotProduct(beam.direction, normals[i]);

		if (cosIncident >= EPS_COS_90) /// beam incidents to this facet
		{
			facetIndices[indicesNumber] = i;
			++indicesNumber;
		}
	}
}

void TracingConcave::SelectVisibleFacetsInternal(const Beam &beam, int *facetIndices,
												 int &indicesNumber)
{
	Point3f *normals = (beam.isExternal) ? m_particle->normals
										 : m_particle->externalNormals;

	Point3f cob = CenterOfPolygon(beam.polygon, beam.size);

	for (int i = 0; i < m_particle->facetNum; ++i)
	{
		double cosIncident = DotProduct(beam.direction, normals[i]);

		if (cosIncident >= EPS_COS_90) /// beam incidents to this facet
		{
			Point3f cof = CenterOfPolygon(m_particle->facets[i], m_particle->vertexNums[i]);
			Point3f vectorToFacet = cof - cob;
			double cosFacets = DotProduct(-normals[beam.facetId]/*OPT:опр.др.норм.*/, vectorToFacet);

			if (cosFacets >= 0.0001/*TODO: подобрать норм значение и вынести*/) /// facet is in front of begin of beam
			{
				facetIndices[indicesNumber] = i;
				++indicesNumber;
			}
		}
	}
}

void TracingConcave::SetPolygonByFacet(const Point3f *facet, int size, Paths &polygon) const
{
	for (int i = 0; i < size; ++i)
	{
		const Point3f &p = facet[i];

		polygon[0] << IntPoint((cInt)std::round(p.cx * MULTI_INDEX),
							   (cInt)std::round(p.cy * MULTI_INDEX),
							   (cInt)std::round(p.cz * MULTI_INDEX));
	}
}

void TracingConcave::CutShadowsFromFacet(const Point3f *facet, int size,
										 int *facetIds, int previewFacetCount,
										 const Beam &beam,
										 Paths &resultPolygon)
{	// TODO: вызывать только если нужно (поставить проверку)
//	int originFacetId = facetIds[previewFacetCount]; /// facet id of current facet
//	Point3f &normal = m_particle->externalNormals[originFacetId];
	Point3f *normals = (beam.isExternal) ? m_particle->externalNormals
										 : m_particle->normals;
	Point3f &originNormal = normals[beam.facetId];

	SetPolygonByFacet(facet, size, resultPolygon); /// set origin polygon

	if (fabs(originNormal.cx) > 0.5)
	{
		ExchangeCoords(Axis::aX, Axis::aZ, resultPolygon);
	}
	else if (fabs(originNormal.cy) > 0.5)
	{
		ExchangeCoords(Axis::aY, Axis::aZ, resultPolygon);
	}

	Paths clip(previewFacetCount);

	for (int i = 0; i < previewFacetCount; ++i)
	{
		// set clip polygon by projection
		int facetId = facetIds[i];
		const Point3f *fac = m_particle->facets[facetId];
		int size = m_particle->vertexNums[facetId];
		ProjectFacetToFacet(fac, size, beam.direction, originNormal, clip[i]);

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

	if (fabs(originNormal.cx) > 0.5)
	{
		ExchangeCoords(Axis::aX, Axis::aZ, clip);
	}
	else if (fabs(originNormal.cy) > 0.5)
	{
		ExchangeCoords(Axis::aY, Axis::aZ, clip);
	}

	Paths result;
	m_clipper.AddPaths(resultPolygon, ptSubject, true);
	m_clipper.AddPaths(clip, ptClip, true);
	m_clipper.Execute(ctDifference, result);
	m_clipper.Clear();

	RemoveEmptyPolygons(result);

	if (result.empty())
	{
		resultPolygon.clear();
		return;
	}
	else /*if (result.size() == MAX_POLYGON_RESULT
			 && result.at(0).size() >= MIN_VERTEX_NUM)*/
	{
		// обратно
		if (fabs(originNormal.cx) > 0.5)
		{
			ExchangeCoords(Axis::aZ, Axis::aX, result);
		}
		else if (fabs(originNormal.cy) > 0.5)
		{
			ExchangeCoords(Axis::aZ, Axis::aY, result);
		}

		ClipperLib::CleanPolygons(result, EPS_MULTI);
		RemoveEmptyPolygons(result);

		resultPolygon = result;
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

void TracingConcave::SetBeamShapeByPolygon(Beam &beam, const Path &result)
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

//	if (beam.size <= 0 || beam.size > MAX_VERTEX_NUM) // DEB
	assert(beam.size > 0 && beam.size <= MAX_VERTEX_NUM);
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
//		Point3f p;
//		ProjectPointToFacet(a_facet[i], a_dir, b_normal, p);

		Point3d p;
		ProjectPointToFacet(Point3d(a_facet[i].cx, a_facet[i].cy, a_facet[i].cz),
							Point3d(a_dir.cx, a_dir.cy, a_dir.cz),
							Point3d(b_normal.cx, b_normal.cy, b_normal.cz, b_normal.d_param), p);

		projection << IntPoint((cInt)std::round(p.x * MULTI_INDEX),
							   (cInt)std::round(p.y * MULTI_INDEX),
							   (cInt)std::round(p.z * MULTI_INDEX));
	}
}

double TracingConcave::MeasureMinDistanceToFacet(int facetId, const Point3f &beamDir)
{
	double dist = FLT_MAX;
	Point3f *facet = m_particle->facets[facetId];

	for (int i = 0; i < m_particle->vertexNums[facetId]; ++i)
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
		distances[i] = MeasureMinDistanceToFacet(facetIds[i], beamDir);
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

void TracingConcave::ExchangeCoords(Axis oldAxis, Axis newAxis, Paths &origin) const
{
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

void TracingConcave::CutBeamShapeByFacet(Paths &beam, int facetId, const Point3f &direction,
										 const Point3f &shapeNormal,
										 Paths &result)
{
	if (fabs(shapeNormal.cx) > 0.5)
	{
		ExchangeCoords(Axis::aX, Axis::aZ, beam);
	}
	else if (fabs(shapeNormal.cy) > 0.5)
	{
		ExchangeCoords(Axis::aY, Axis::aZ, beam);
	}

	Paths clip(1);
	{
		const Point3f *facet = m_particle->facets[facetId];
		int size = m_particle->vertexNums[facetId];
		ProjectFacetToFacet(facet, size, direction, shapeNormal, clip[0]); // проецируем грань на начальный пучок
	}

	if (fabs(shapeNormal.cx) > 0.5)
	{
		ExchangeCoords(Axis::aX, Axis::aZ, clip);
	}
	else if (fabs(shapeNormal.cy) > 0.5)
	{
		ExchangeCoords(Axis::aY, Axis::aZ, clip);
	}

	m_clipper.AddPaths(beam, ptSubject, true);
	m_clipper.AddPaths(clip, ptClip, true);
	m_clipper.Execute(ctDifference, result);
	m_clipper.Clear();

	if (!result.empty())
	{
		// обратно
		if (fabs(shapeNormal.cx) > 0.5)
		{
			ExchangeCoords(Axis::aZ, Axis::aX, result);
		}
		else if (fabs(shapeNormal.cy) > 0.5)
		{
			ExchangeCoords(Axis::aZ, Axis::aY, result);
		}

		ClipperLib::CleanPolygons(result, EPS_MULTI);
		RemoveEmptyPolygons(result);

//		if (result.size() >= 2)
//		{
//			int fff = 0;
//			ClipperLib::CleanPolygon(result[0], EPS_MULTI);
//		}
	}
}

void TracingConcave::RemoveEmptyPolygons(Paths &result)
{
	Paths buff = result; // OPT: переделать
	result.clear();

	for (int i = 0; i < buff.size(); ++i)
	{
		if (buff.at(i).size() >= MAX_POLYGON_RESULT)
		{
			result.push_back(buff.at(i));
		}
	}
}

void TracingConcave::CutBeamShapeByFacet(int facetId, const Beam &beam,
										 const Point3f &shapeNormal,
										 Paths &result)
{
	Paths origin(1);
	SetPolygonByFacet(beam.polygon, beam.size, origin);

	if (fabs(shapeNormal.cx) > 0.5)
	{
		ExchangeCoords(Axis::aX, Axis::aZ, origin);
	}
	else if (fabs(shapeNormal.cy) > 0.5)
	{
		ExchangeCoords(Axis::aY, Axis::aZ, origin);
	}

	Paths clip(1);
	{
		const Point3f *facet = m_particle->facets[facetId];
		int size = m_particle->vertexNums[facetId];
		ProjectFacetToFacet(facet, size, beam.direction, shapeNormal, clip[0]); // проецируем грань на начальный пучок
	}

	if (fabs(shapeNormal.cx) > 0.5)
	{
		ExchangeCoords(Axis::aX, Axis::aZ, clip);
	}
	else if (fabs(shapeNormal.cy) > 0.5)
	{
		ExchangeCoords(Axis::aY, Axis::aZ, clip);
	}

	m_clipper.AddPaths(origin, ptSubject, true);
	m_clipper.AddPaths(clip, ptClip, true);
	m_clipper.Execute(ctDifference, result);
	m_clipper.Clear();

	if (!result.empty())
	{
		// обратно
		if (fabs(shapeNormal.cx) > 0.5)
		{
			ExchangeCoords(Axis::aZ, Axis::aX, result);
		}
		else if (fabs(shapeNormal.cy) > 0.5)
		{
			ExchangeCoords(Axis::aZ, Axis::aY, result);
		}

		ClipperLib::CleanPolygons(result, EPS_MULTI);
		RemoveEmptyPolygons(result);
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

	assert(false && "Divide point is not found");
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
		Point3f v1 = polygon[i-1] - polygon[baseI]; // OPT:
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

	if (square < 0)
	{	/// OPT: узнать в какую сторону ориентированы точки в пучке
		square *= (-1);
	}

	return square/2;
}

void MeasureZ(IntPoint &a1, IntPoint &a2, IntPoint &, IntPoint &, IntPoint &point)
{
	IntPoint top, bot; // OPT: реализовать на указателях

	if (a2.Z > a1.Z)
	{
		top = a2;
		bot = a1;
	}
	else
	{
		top = a1;
		bot = a2;
	}

	double x, y;

	IntPoint vec;
	vec.X = top.X - bot.X;
	vec.Y = top.Y - bot.Y;

	x = vec.X;
	y = vec.Y;
	double normVec = x*x + y*y;

	IntPoint botVec;
	botVec.X = point.X - bot.X;
	botVec.Y = point.Y - bot.Y;

	x = botVec.X;
	y = botVec.Y;
	double normBotVec = x*x + y*y;

	point.Z = (top.Z - bot.Z) * sqrt(normBotVec) / sqrt(normVec) + bot.Z;
	int fff = 0;
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
		ExchangeCoords(Axis::aX, Axis::aZ, polygon);
		float tmp = normal1.cx;
		normal1.cx = normal1.cz;
		normal1.cz = tmp;
	}
	else if (fabs(n_normal.cy) > 0.5)
	{
		ExchangeCoords(Axis::aY, Axis::aZ, polygon);
		float tmp = normal1.cy;
		normal1.cy = normal1.cz;
		normal1.cz = tmp;
	}

	area = ClipperLib::Area(polygon[0]);
	area /= (cInt)MULTI_INDEX * MULTI_INDEX; // OPT

	Point3f z(0, 0, 1);
	double cosNormZ = DotProduct(normal1, z);

	if (cosNormZ > 0)
	{
		cosNormZ = DotProduct(normal1, -z); // OPT ?
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
	const double eps = 1e7*DBL_EPSILON; // OPT

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
	double n = sqrt(Norm(normal));
	return (e*square) / n; // REF сделать функ length
}
