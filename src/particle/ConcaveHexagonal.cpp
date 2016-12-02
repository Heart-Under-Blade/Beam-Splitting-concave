#include "ConcaveHexagonal.h"

ConcaveHexagonal::ConcaveHexagonal(double cavityDept)
{
	m_cavityDept = cavityDept;
}

ConcaveHexagonal::ConcaveHexagonal(double radius, double halfHeight, const complex &refractionIndex,
								   double cavityDept)
	: Hexagonal(radius, halfHeight, refractionIndex)
{
	m_cavityDept = cavityDept;

	SetFacetParams();

	SetOriginBaseFacets();
	SetOriginCavityPoints();

	SetCavities(m_originBases.top, m_originBases.bottom, m_originCavities);
	SetSideFacets(m_originBases.top, m_originBases.bottom,
				  BASE_VERTEX_NUM, 2*BASE_VERTEX_NUM);

	SetOriginNormals();
}

void ConcaveHexagonal::Rotate(double beta, double gamma, double alpha)
{
	SetRotateMatrix(beta, gamma, alpha);

	Point3f baseTop[BASE_VERTEX_NUM], baseBottom[BASE_VERTEX_NUM];
	RotateBaseFacets(baseTop, baseBottom);

	CavityPoints cavities;
	RotatePoint(m_originCavities.top, cavities.top);
	RotatePoint(m_originCavities.bottom, cavities.bottom);

	SetCavities(baseTop, baseBottom, cavities);
	SetSideFacets(baseTop, baseBottom, BASE_VERTEX_NUM, 2*BASE_VERTEX_NUM);

	RotateNormals();
}

void ConcaveHexagonal::SetFacetParams()
{
	facetNum = BASE_VERTEX_NUM + BASE_FACET_NUM * BASE_VERTEX_NUM;

	// top facet (triangles)
	for (int i = 0; i < BASE_VERTEX_NUM; ++i)
	{
		vertexNums[i] = CAVITY_FACET_VERTEX_NUM;
	}

	// side facet
	for (int i = BASE_VERTEX_NUM; i < 2*BASE_VERTEX_NUM; ++i)
	{
		vertexNums[i] = SIDE_VERTEX_NUMBER;
	}

	// bottom facet (triangles)
	for (int i = 2*BASE_VERTEX_NUM; i < 3*BASE_VERTEX_NUM; ++i)
	{
		vertexNums[i] = CAVITY_FACET_VERTEX_NUM;
	}
}

void ConcaveHexagonal::SetCavities(Point3f *baseTop, Point3f *baseBottom,
								   const CavityPoints &cavities)
{
	SetCavityFacets(0, BASE_VERTEX_NUM, baseTop, cavities.top); // top facets (triangles)
	SetCavityFacets(2*BASE_VERTEX_NUM, 3*BASE_VERTEX_NUM, baseBottom, cavities.bottom); // bottom facets (triangles)
}

void ConcaveHexagonal::SetCavityFacets(int start, int end, Point3f *baseFacet,
									   const Point3f &cavityPoint)
{
	// base facet point indices
	int p0 = BASE_VERTEX_NUM-1;
	int p1 = 0;

	for (int i = start; i < end; ++i)
	{
		facets[i][0] = baseFacet[p0];
		facets[i][1] = baseFacet[p1];
		facets[i][2] = cavityPoint;
		p0 = p1;
		++p1;
	}
}

void ConcaveHexagonal::SetOriginCavityPoints()
{
	m_originCavities.top = Point3f(0, 0, halfHeight - m_cavityDept);
	m_originCavities.bottom = Point3f(0, 0, -halfHeight + m_cavityDept);
}

void ConcaveHexagonal::SetBaseNormals()
{
	// top cavity facets
	for (int i = 0; i < BASE_VERTEX_NUM; ++i)
	{
		m_originNormals[i] = -NormalToFacet(facets[i]);
	}

	// bottom cavity facets
	for (int i = 2*BASE_VERTEX_NUM; i < 3*BASE_VERTEX_NUM; ++i)
	{
		m_originNormals[i] = -NormalToFacet(facets[i]);
	}
}

void ConcaveHexagonal::SetSideNormals()
{
	// side facets
	double cos30 = sqrt(3)/2;
	m_originNormals[6]	= Point3f(-cos30, -0.5, 0); // REF: то же самое есть в базовом классе
	m_originNormals[7]	= Point3f(0, -1, 0);
	m_originNormals[8]	= Point3f(cos30, -0.5, 0);
	m_originNormals[9]	= Point3f(cos30, 0.5, 0);
	m_originNormals[10] = Point3f(0, 1, 0);
	m_originNormals[11] = Point3f(-cos30, 0.5, 0);
}

void ConcaveHexagonal::SetOriginNormals()
{
	SetBaseNormals();
	SetSideNormals();

	CopyPoints(m_originNormals, normals, facetNum); // current normals

	SetDParams();
	SetExternalNormals();
}
