#pragma once

#include "Tracing.h"
#include "clipper.hpp"

typedef std::vector<std::vector<Point3f>> Polygons;

struct BeamInfoConcave : public BeamInfo
{// OPT: попробовать объеденить с Beam
	bool isExternal;

	BeamInfoConcave() {}
	BeamInfoConcave(Beam b, int fi, int d, bool ie)
	{
		beam = b;
		facetId = fi;
		dept = d;
		isExternal = ie;
	}
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

private:
	enum class Axis : int
	{
		aX, aY, aZ
	};

	ClipperLib::Clipper m_clipper;

private:
	double MeasureMinDistanceToFacet(int facetId, const Point3f &beamDir);
	void SortFacets(int number, const Point3f &beamDir, int *facetIds); ///< use fast sort algorithm
	void SelectVisibleFacetsExternal(const Beam &beam, int *facetIndices,
									 int &indicesNumber);
	void SelectVisibleFacetsInternal(const BeamInfoConcave &beamInfo, int *facetIndices,
									 int &indicesNumber);
	void CutFacetByFacetShadows(int *facetIds, const Point3f *facet, int size, int previewFacetCount,
						   ClipperLib::Paths &resultFacet);

	void ProjectPointToFacet(const Point3f &point, const Point3f &direction,
							 const Point3f &facetNormal, Point3f &projection);
	void ProjectFacetToFacet(const Point3f *a_facet, int a_size, const Point3f &b_normal,
							 ClipperLib::Paths &projection);

	void SetBeamShapeByPolygon(Beam &beam, const ClipperLib::Paths &result);

	void CutBeamShapeByFacet(int facetId, const Beam &beam, const Point3f &shapeNormal,
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

	void ExchangeCoords(Axis oldAxis, Axis newAxis, ClipperLib::Paths &origin); ///< заменяем координаты, для устранения погрешности при клиппинге

	void SetPolygonByFacet(const Point3f *facet, int size, ClipperLib::Paths &polygon) const;

	void CutReflectedBeam(const BeamInfoConcave &beamInfo, Beam &incidentBeam);
        
protected:
	void TraceInternalReflections(BeamInfoConcave *tree, int size,
								  std::vector<Beam> &outBeams);
};

void MeasureZ(ClipperLib::IntPoint & a1, ClipperLib::IntPoint & a2,
			  ClipperLib::IntPoint &, ClipperLib::IntPoint &,
			  ClipperLib::IntPoint & point);
