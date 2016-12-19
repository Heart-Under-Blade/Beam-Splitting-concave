#include "Particle.h"

Particle::Particle() {}

Particle::Particle(double p_radius, double p_halfHeight,
				   const complex &p_refractionIndex)
{
	Init(p_radius, p_halfHeight, p_refractionIndex);
}

double Particle::GetHalfHeight() const
{
    return m_halfHeight;
}

void Particle::Init(double p_radius, double p_halfHeight, const complex &p_refractionIndex)
{
	m_radius = p_radius;
	m_halfHeight = p_halfHeight;
	refractionIndex = p_refractionIndex;

	double re = real(refractionIndex);
	double im = imag(refractionIndex);
	ri_coef_re = re*re - im*im;
	ri_coef_im = 4*re*re*im;
}

void Particle::SetRotateMatrix(double beta, double gamma, double alpha)
{
	double cosA, cosB, cosG,
			sinA, sinB, sinG;

	sincos(alpha, &sinA, &cosA);
	sincos(beta, &sinB, &cosB);
	sincos(gamma, &sinG, &cosG);

	double cosAcosB = cosA*cosB;
	double sinAcosG = sinA*cosG;
	double sinAsinG = sinA*sinG;

	m_rotMatrix[0][0] = cosAcosB*cosG - sinAsinG;
	m_rotMatrix[1][0] = sinAcosG*cosB + cosA*sinG;
	m_rotMatrix[2][0] = -sinB*cosG;

	m_rotMatrix[0][1] = -(cosAcosB*sinG + sinAcosG);
	m_rotMatrix[1][1] = cosA*cosG - sinAsinG*cosB;
	m_rotMatrix[2][1] = sinB*sinG;

	m_rotMatrix[0][2] = cosA*sinB;
	m_rotMatrix[1][2] = sinA*sinB;
	m_rotMatrix[2][2] = cosB;
}

void Particle::RotateNormals()
{
	for (int i = 0; i < facetNum; ++i)
	{
		RotatePoint(m_originNormals[i], facets[i].in_normal);
	}

	SetDParams();
	SetExternalNormals();
}

void Particle::SetDParams()
{
	for (int i = 0; i < facetNum; ++i)
	{
		double d = DotProduct(facets[i].polygon.arr[0], facets[i].in_normal);
		facets[i].in_normal.d_param = -d;
	}
}

void Particle::SetExternalNormals()
{
	for (int i = 0; i < facetNum; ++i)
	{
		Facet &facet = facets[i];
		facet.ex_normal.cx = -facet.in_normal.cx;
		facet.ex_normal.cy = -facet.in_normal.cy;
		facet.ex_normal.cz = -facet.in_normal.cz;
		facet.ex_normal.d_param = -facet.in_normal.d_param;
	}
}

void Particle::CopyFacet(Point3f *points, Facet &result)
{
	for (int i = 0; i <= result.polygon.size; ++i)
	{
		result.polygon.arr[i].cx = points[i].cx;
		result.polygon.arr[i].cy = points[i].cy;
		result.polygon.arr[i].cz = points[i].cz;
	}
}

void Particle::RotatePoint(const Point3f &point, Point3f &result)
{
	result.cx = point.cx*m_rotMatrix[0][0] + point.cy*m_rotMatrix[0][1] + point.cz*m_rotMatrix[0][2];
	result.cy = point.cx*m_rotMatrix[1][0] + point.cy*m_rotMatrix[1][1] + point.cz*m_rotMatrix[1][2];
	result.cz = point.cx*m_rotMatrix[2][0] + point.cy*m_rotMatrix[2][1] + point.cz*m_rotMatrix[2][2];
}
