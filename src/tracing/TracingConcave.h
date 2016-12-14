#pragma once

#include "Tracing.h"
#include "clipper.hpp"

typedef std::vector<std::vector<Point3f>> Polygons;

enum class Axis : int
{
	aX, aY, aZ
};

/** NOTE: пучки выходят со случайно ориентированными вершинами */
class TracingConcave : public Tracing
{
public:
	TracingConcave(Particle *particle, const Point3f &startBeamDir,
				   bool isOpticalPath, const Point3f &polarizationBasis,
				   int interReflectionNumber);

	double BeamCrossSection(const Beam &beam) const override;

	void SplitBeamByParticle(std::vector<Beam> &outBeams,
							 double &lightSurfaceSquare) override;

	void SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
							 std::vector<Beam> &outBeams) override;

	void CutIncidentBeam(int facetId, Beam &beam, bool &isDivided);

private:
	ClipperLib::Clipper m_clipper;

private:

	double CalcMinDistanceToFacet(int facetId, const Point3f &beamDir);
	void SortFacets(int number, const Point3f &beamDir, int *facetIds); ///< use fast sort algorithm
	void CutShadowsFromFacet(const Point3f *facet, int size, int *facetIds,
								int previewFacetCount, const Beam &incidentBeam,
								ClipperLib::Paths &resultFacet);

	void ProjectPointToFacet(const Point3d &point, const Point3d &direction,
							 const Point3d &facetNormal, Point3d &projection);

	void ProjectPointToFacet(const Point3f &point, const Point3f &direction,
							 const Point3f &facetNormal, Point3f &projection);

	void ProjectFacetToFacet(const Point3f *a_facet, int a_size, const Point3f &a_dir, const Point3f &b_normal,
							 ClipperLib::Path &projection);

	void SetBeamByPath(Beam &beam, const ClipperLib::Path &result);

	void CutBeamByFacet(ClipperLib::Paths &beamPolygon, int facetId,
							 const Point3f &direction, const Point3f &polygonNormal,
							 ClipperLib::Paths &result);

	void DivideConcavePolygon(const Point3f *polygon, int size,
							  const Point3f &normal,
							  Polygons &polygons) const;

	void FindDividePoint(const std::vector<Point3f> &polygon,
						 int i0, int i1, const Point3f &normal,
						 Point3f &x, int &nextPointIndex) const;

	void FillSubpolygon(int begin, int end,
						const std::vector<Point3f> &polygon,
						std::vector<Point3f> &subpolygon) const;

	// recursive
	void DividePolygon(const std::vector<Point3f> &polygon,
					   const Point3f &normal, Polygons &polygons) const;

	double SquareOfPolygon(const std::vector<Point3f> &polygon) const;

	void SwapCoords(Axis oldAxis, Axis newAxis, ClipperLib::Paths &origin) const; ///< заменяем координаты, для устранения погрешности при клиппинге

	void SetPolygonByFacet(const Point3f *facet, int size, ClipperLib::Paths &polygon) const;
        
	double AreaByClipper(const Beam &beam, const Point3f &normal) const;

//	bool isOrderReversed(const Point3f oldNormal, const ClipperLib::Path polygon);

	void InversePolygonOrder(ClipperLib::Path &polygon);

	void CatchExternalBeam(const Beam &beam, std::vector<Beam> &outBeams);

	void PushBeamToTree(Beam &beam, int facetId, int level, bool isExternal);

	void SelectVisibleFacetsExternal(const Beam &beam, int *facetIndices,
									 int &indicesNumber);
	void FindVisibleFacetsInternal(const Beam &beam, int *facetIndices,
								   int &facetIdCount);
	void RemoveEmptyPolygons(ClipperLib::Paths &result);
        
	void PrintTrack(const Beam &beam, int facetId);

	Axis GetSwapAxis(const Point3f &normal);
        
	void ClipDifference(const ClipperLib::Paths &subject, const ClipperLib::Paths &clip,
						ClipperLib::Paths &difference);

	void SelectVisibleFacets(const Beam &beam, int *facetIds, int &facetIdCount);

	void RemoveHole(ClipperLib::Paths &result);

	void HandleResultPolygon(Axis axis, ClipperLib::Paths &result);

protected:
	void TraceInternalReflections(std::vector<Beam> &outBeams);
};

void FindZCoord(ClipperLib::IntPoint & a1, ClipperLib::IntPoint & a2,
			  ClipperLib::IntPoint &, ClipperLib::IntPoint &,
			  ClipperLib::IntPoint & point);
