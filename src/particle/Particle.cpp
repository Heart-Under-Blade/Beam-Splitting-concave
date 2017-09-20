#include "Particle.h"
#include <fstream>

Particle::Particle() {}

void Particle::SetFromFile(const char *filename)
{
	std::ifstream pfile(filename, std::ios::in);
	//pfile >>
}

void Particle::Init(int facetCount, const complex &refrIndex, double size)
{
	facetNum = facetCount;
	m_refractiveIndex = refrIndex;
	m_mainSize = size;
}

void Particle::RotateCenters()
{
	for (int i = 0; i < facetNum; ++i)
	{
		RotatePoint(defaultFacets[i].center, facets[i].center);
	}
}

void Particle::Rotate(double beta, double gamma, double alpha)
{
	SetRotateMatrix(beta, gamma, alpha);

	// REF: слить всё в один цикл
	for (int i = 0; i < facetNum; ++i)
	{
		for (int j = 0; j < facets[i].size; ++j)
		{
			RotatePoint(defaultFacets[i].arr[j], facets[i].arr[j]);
		}
	}

	RotateNormals();
	RotateCenters();
}

void Particle::Fix()
{
	for (int i = 0; i < facetNum; ++i)
	{
		for (int j = 0; j < facets[i].size; ++j)
		{
			defaultFacets[i].arr[j] = facets[i].arr[j];
		}
	}
}

void Particle::Concate(const std::vector<Particle> &parts)
{
	int i = 0;
	facetNum = 0;

	for (const Particle &part : parts)
	{
		facetNum += part.facetNum;

		for (int j = 0; j < part.facetNum; ++j)
		{
			defaultFacets[i++] = part.facets[j];
		}
	}

	isAggregate = true;
}

const double &Particle::GetMainSize() const
{
	return m_mainSize;
}

const complex &Particle::GetRefractionIndex() const
{
	return m_refractiveIndex;
}

const Symmetry &Particle::GetSymmetry() const
{
	return m_symmetry;
}

void Particle::Move(float dx, float dy, float dz)
{
	for (int i = 0; i < facetNum; ++i)
	{
		for (int j = 0; j < defaultFacets[i].size; ++j)
		{
			facets[i].arr[j] = defaultFacets[i].arr[j] + Point3f(dx, dy, dz);
		}
	}
}

void Particle::Output()
{
	std::ofstream M("particle.dat", std::ios::out);

	for (int i = 0; i < facetNum; ++i)
	{
		for (int j = 0; j < facets[i].size; ++j)
		{
			Point3f p = facets[i].arr[j];
			M << p.point[0] << ' '
							<< p.point[1] << ' '
							<< p.point[2] << ' '
							<< i ;
			M << std::endl;
		}

		M << std::endl << std::endl;;
	}

	M.close();
}


void Particle::SetDefaultNormals()
{
	for (int i = 0; i < facetNum; ++i)
	{
		defaultFacets[i].SetNormal();
	}
}

void Particle::Reset()
{
	for (int i = 0; i < facetNum; ++i)
	{
		facets[i] = defaultFacets[i];
	}
}

void Particle::SetDefaultCenters()
{
	for (int i = 0; i < facetNum; ++i)
	{
		defaultFacets[i].SetCenter();
	}
}

void Particle::SetRotateMatrix(double beta, double gamma, double alpha)
{
	double cosA, cosB, cosG,
			sinA, sinB, sinG;

	sincos(alpha, &sinA, &cosA);
	sincos(beta,  &sinB, &cosB);
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
		RotatePoint(defaultFacets[i].in_normal, facets[i].in_normal);
	}

	SetDParams();

	for (int i = 0; i < facetNum; ++i)
	{
		facets[i].ex_normal = -facets[i].in_normal;
		facets[i].ex_normal.d_param = -facets[i].in_normal.d_param;
	}
}

void Particle::SetDParams()
{
	for (int i = 0; i < facetNum; ++i)
	{
		double d = DotProduct(facets[i].arr[0], facets[i].in_normal);
		facets[i].in_normal.d_param = -d;
	}
}

void Particle::RotatePoint(const Point3f &point, Point3f &result)
{
	result.cx = point.cx*m_rotMatrix[0][0] + point.cy*m_rotMatrix[0][1] + point.cz*m_rotMatrix[0][2];
	result.cy = point.cx*m_rotMatrix[1][0] + point.cy*m_rotMatrix[1][1] + point.cz*m_rotMatrix[1][2];
	result.cz = point.cx*m_rotMatrix[2][0] + point.cy*m_rotMatrix[2][1] + point.cz*m_rotMatrix[2][2];
}

void Particle::SetSymmetry(double beta, double gamma, double alpha)
{
	m_symmetry.beta = beta;
	m_symmetry.gamma = gamma;
	m_symmetry.alpha = alpha;
}
