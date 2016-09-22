#include "Particle.h"

Particle::Particle(double p_radius, double p_halfHeight,
				   const complex &p_refractionIndex)
{
	m_radius = p_radius;
	m_halfHeight = p_halfHeight;
	refractionIndex = p_refractionIndex;

	double re = real(refractionIndex);
	double im = imag(refractionIndex);
	double re_sqr = re*re;
	refrI_sqr_re = re_sqr - im*im;
	refrI_sqr_im = 4*re_sqr*im;
}

void Particle::SetRotateMatrix(double beta, double gamma, double alpha)
{
	double cosA, cosB, cosG,
			sinA, sinB, sinG;

#ifdef __FASTMATH_H
	sincos(alpha, &sinA, &cosA);
	sincos(beta, &sinB, &cosB);
	sincos(gamma, &sinG, &cosG);
#else
	// TODO: опт. заменить на sincos'ы
	cosA = cos(alpha);
	cosB = cos(beta);
	cosG = cos(gamma);

	sinA = sin(alpha);
	sinB = sin(beta);
	sinG = sin(gamma);
#endif

	/// TODO: дооптимизировать
	double cosAcosB = cosA*cosB;

	m_rotMatrix[0][0] = cosAcosB*cosG - sinA*sinG;
	m_rotMatrix[1][0] = sinA*cosB*cosG + cosA*sinG;
	m_rotMatrix[2][0] = -sinB*cosG;

	m_rotMatrix[0][1] = -(cosAcosB*sinG + sinA*cosG);
	m_rotMatrix[1][1] = cosA*cosG - sinA*cosB*sinG;
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
		normals[i].D_PARAM = -d;
	}
}

void Particle::SetExternalNormals()
{
	for (int i = 0; i < facetNum; ++i)
	{
		externalNormals[i].X = -normals[i].X;
		externalNormals[i].Y = -normals[i].Y;
		externalNormals[i].Z = -normals[i].Z;
		externalNormals[i].D_PARAM = -normals[i].D_PARAM;
	}
}

void Particle::RotatePoint(const Point3f &point, Point3f &result)
{
	result.X = point.X*m_rotMatrix[0][0] + point.Y*m_rotMatrix[0][1] + point.Z*m_rotMatrix[0][2];
	result.Y = point.X*m_rotMatrix[1][0] + point.Y*m_rotMatrix[1][1] + point.Z*m_rotMatrix[1][2];
	result.Z = point.X*m_rotMatrix[2][0] + point.Y*m_rotMatrix[2][1] + point.Z*m_rotMatrix[2][2];
}
