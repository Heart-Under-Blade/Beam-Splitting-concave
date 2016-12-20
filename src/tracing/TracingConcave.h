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

	void SplitBeamByParticle(std::vector<Beam> &outBeams) override;

	void SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
							 std::vector<Beam> &outBeams) override;

	void CutIncidentBeam(int facetId, Beam &beam, bool &isDivided);

private:
	ClipperLib::Clipper m_clipper;

private:
	double CalcMinDistanceToFacet(const Polygon &polygon, const Point3f &beamDir);
	void SortFacets(const Point3f &beamDir, IntArray &facetIds); ///< use fast sort algorithm
	void CutShadowsFromFacet(int facetId, IntArray facetIds, int previewFacetCount,
							 const Beam &incidentBeam, ClipperLib::Paths &resultFacet);

	void ProjectPointToFacet(const Point3d &point, const Point3d &direction,
							 const Point3d &facetNormal, Point3d &projection);

	void ProjectPointToFacet(const Point3f &point, const Point3f &direction,
							 const Point3f &facetNormal, Point3f &projection);

	void ProjectFacetToFacet(const Point3f *a_facet, int a_size, const Point3f &a_dir, const Point3f &b_normal,
							 ClipperLib::Path &projection);

	void PathToPolygon(const ClipperLib::Path &path, Polygon &polygon);

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

	void SwapCoords(Axis oldAxis, Axis newAxis, ClipperLib::Paths &origin) const; ///< заменяем координаты, для устранения погрешности при клиппинге

	void PolygonToPath(const Polygon &pol, ClipperLib::Paths &path) const;
        
	double AreaByClipper(const Beam &beam, const Point3f &normal) const;

//	bool isOrderReversed(const Point3f oldNormal, const ClipperLib::Path polygon);

//	void InversePolygonOrder(ClipperLib::Path &polygon);

	void CatchExternalBeam(const Beam &beam, std::vector<Beam> &scatteredBeams);

	void PushBeamToTree(Beam &beam, int facetId, int level, Location location);

	void FindVisibleFacets_initial(const Beam &beam, IntArray &facetIds);
	void FindVisibleFacets(const Beam &beam, IntArray &facetIds);
	void RemoveEmptyPolygons(ClipperLib::Paths &result);

	Axis GetSwapAxis(const Point3f &normal);
        
	void ClipDifference(const ClipperLib::Paths &subject, const ClipperLib::Paths &clip,
						ClipperLib::Paths &difference);

	void SelectVisibleFacets(const Beam &beam, IntArray &facetIds);

	void RemoveHole(ClipperLib::Paths &result);

	void HandleResultPolygon(Axis axis, ClipperLib::Paths &result);

	void SetOpticalBeamParams(int facetId, Beam &incidentBeam,
							  Beam &inBeam, Beam &outBeam, bool &hasOutBeam);

	void SetIntersectedPolygon(const IntArray &facetIds, int handledFacetNum,
							   Polygon &beamPolygon, bool &isTotallyShadowed);

	void TraceInitialBeam();

	void SelectVisibleFacets_initial(IntArray &facetIds);

	bool IsRemainedExternalBeam(Beam &incidentBeam);

#ifdef _TRACK_ALLOW
	void PrintTrack(const Beam &beam, int facetId);
	void AddToTrack(Beam &beam, int facetId);
#endif

protected:
	void TraceSecondaryBeams(std::vector<Beam> &scaterredBeams);
};

void FindZCoord(ClipperLib::IntPoint & a1, ClipperLib::IntPoint & a2,
			  ClipperLib::IntPoint &, ClipperLib::IntPoint &,
			  ClipperLib::IntPoint & point);
