#pragma once

#include "clipper.hpp"
#include "geometry_lib.h"

#define CLIP_RESULT_SINGLE 1

#define MULTI_INDEX		10000000l	// index for poligon's clip operations
#define EPS_MULTI		(1.415*MULTI_INDEX*2)/10000	// погрешность, при которой точки операций Clipper'а можно считать совпадающими

enum class Axis : int
{
	aX, aY, aZ
};

class BeamClipper
{
public:
	BeamClipper();

	Axis GetSwapAxis(const Point3f &normal) const;

	/**
	 * @brief SwapCoords Заменяем координаты, для устранения погрешности при клиппинге
	 * @param oldAxis Первая ось обмена координат
	 * @param newAxis Вторая ось обмена координат
	 * @param origin Результирующий полигон
	 */
	void SwapCoords(Axis ax1, Axis ax2, ClipperLib::Paths &origin) const;

	void PolygonToPath(const Polygon &pol, ClipperLib::Paths &path) const;
	void PathToPolygon(const ClipperLib::Path &path, Polygon &polygon) const;

	double AreaOfConcavePolygon(const Polygon &beam, const Point3f &normal) const;

	void Difference(const ClipperLib::Paths &subject, const ClipperLib::Paths &clip,
						ClipperLib::Paths &difference);

	void RemoveHole(ClipperLib::Paths &result) const;
	void RemoveEmptyPaths(ClipperLib::Paths &result) const;
	void HandleResultPaths(Axis axis, ClipperLib::Paths &result) const;

private:
	ClipperLib::Clipper m_clipper;
};

void FindZCoord(ClipperLib::IntPoint & a1, ClipperLib::IntPoint & a2,
				ClipperLib::IntPoint &, ClipperLib::IntPoint &,
				ClipperLib::IntPoint & point);
