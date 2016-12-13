#include "Particle.h"

Particle::Particle() {}

Particle::Particle(double p_radius, double p_halfHeight,
				   const complex &p_refractionIndex)
{
	Init(p_radius, p_halfHeight, p_refractionIndex);
}

void Particle::Init(double p_radius, double p_halfHeight, const complex &p_refractionIndex)
{
	m_radius = p_radius;
	halfHeight = p_halfHeight;
	refractionIndex = p_refractionIndex;

	double re = real(refractionIndex);
	double im = imag(refractionIndex);
	refrI_coef_re = re*re - im*im;
	refrI_coef_im = 4*re*re*im;
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
		RotatePoint(m_originNormals[i], normals[i]);
	}

	SetDParams();
	SetExternalNormals();
}

void Particle::SetDParams()
{
	for (int i = 0; i < facetNum; ++i)
	{
		double d = DotProduct(facets[i][0], normals[i]);
		normals[i].d_param = -d;
	}
}

void Particle::SetExternalNormals()
{
	for (int i = 0; i < facetNum; ++i)
	{
		externalNormals[i].cx = -normals[i].cx;
		externalNormals[i].cy = -normals[i].cy;
		externalNormals[i].cz = -normals[i].cz;
		externalNormals[i].d_param = -normals[i].d_param;
	}
}

void Particle::RotatePoint(const Point3f &point, Point3f &result)
{
	result.cx = point.cx*m_rotMatrix[0][0] + point.cy*m_rotMatrix[0][1] + point.cz*m_rotMatrix[0][2];
	result.cy = point.cx*m_rotMatrix[1][0] + point.cy*m_rotMatrix[1][1] + point.cz*m_rotMatrix[1][2];
	result.cz = point.cx*m_rotMatrix[2][0] + point.cy*m_rotMatrix[2][1] + point.cz*m_rotMatrix[2][2];
}
