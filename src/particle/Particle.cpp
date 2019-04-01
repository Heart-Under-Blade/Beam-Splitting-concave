#include "Particle.h"
#include <fstream>
#include <iostream>
#include <cstring>
#include "global.h"

Particle::Particle()
{
	isConcave = false;
}

void Particle::SetFromFile(const std::string &filename)
{
	std::ifstream pfile(filename, std::ios::in);

	if (!pfile.is_open())
	{
		std::cerr << "File \"" << filename << "\" is not found" << std::endl;
		throw std::exception();
	}

	const int bufSize = 1024;
	char *buff = (char*)malloc(sizeof(char) * bufSize);

	nFacets = 0;
	Facet *facet = &(defaultFacets[nFacets++]);

	char *ptr, *trash;

	pfile.getline(buff, bufSize);
	isConcave = strtol(buff, &trash, 10);

	pfile.getline(buff, bufSize);
	isAggregated = strtol(buff, &trash, 10);

	// read symmetry params
	{
		pfile.getline(buff, bufSize);

		ptr = strtok(buff, " ");
		m_symmetry.beta = DegToRad(strtod(ptr, &trash));

		ptr = strtok(NULL, " ");
		m_symmetry.gamma = DegToRad(strtod(ptr, &trash));
	}

	pfile.getline(buff, bufSize); // skip empty line

	while (!pfile.eof()) // read vertices of facets
	{
		pfile.getline(buff, bufSize);
		ptr = strtok(buff, " ");

		if (strlen(buff) == 0)
		{
			facet = &(defaultFacets[nFacets++]);
			continue;
		}

		int c_i = 0;

		while (ptr != NULL)
		{
			facet->arr[facet->size].point[c_i++] = strtod(ptr, &trash);
			ptr = strtok(NULL, " ");
		}

		++(facet->size);
	}

	pfile.close();

	// correction of number of facet
	if (defaultFacets[nFacets-1].size == 0)
	{
		--nFacets;
	}

	SetDefaultNormals();
	Reset();
	SetDefaultCenters();

	if (isConcave || isAggregated)
	{
		for (int i = 0; i < nFacets; ++i)
		{
			defaultFacets[i].isVisibleIn = false;
			defaultFacets[i].isVisibleOut = false;
			facets[i].isVisibleIn = false;
			facets[i].isVisibleOut = false;
		}
	}
}

void Particle::Init(int facetCount, const complex &refrIndex)
{
	nFacets = facetCount;
	m_refractiveIndex = refrIndex;
}

void Particle::RotateCenters()
{
	for (int i = 0; i < nFacets; ++i)
	{
		RotatePoint(defaultFacets[i].center, facets[i].center);
	}
}

void Particle::Rotate(double beta, double gamma, double alpha)
{
	rotAngle = Angle{alpha, beta, gamma};
	SetRotateMatrix(beta, gamma, alpha);

	// REF: слить всё в один цикл
	for (int i = 0; i < nFacets; ++i)
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
	for (int i = 0; i < nFacets; ++i)
	{
		for (int j = 0; j < facets[i].size; ++j)
		{
			defaultFacets[i].arr[j] = facets[i].arr[j];
		}
	}
}

void Particle::Scale(double ratio)
{
	for (int i = 0; i < nFacets; ++i)
	{
		for (int j = 0; j < defaultFacets[i].size; ++j)
		{
			defaultFacets[i].arr[j].cx *= ratio;
			defaultFacets[i].arr[j].cy *= ratio;
			defaultFacets[i].arr[j].cz *= ratio;
		}
	}
}

void Particle::Resize(double size)
{
	double radius = GetRotationRadius();
	double ratio = size/radius*2;
	Scale(ratio);
	Reset();
}

void Particle::Concate(const std::vector<Particle> &parts)
{
	int i = 0;
	nFacets = 0;

	for (const Particle &part : parts)
	{
		nFacets += part.nFacets;

		for (int j = 0; j < part.nFacets; ++j)
		{
			defaultFacets[i++] = part.facets[j];
		}
	}

	isAggregated = true;
}

double Particle::GetRotationRadius() const
{
	Point3f p0(0, 0, 0);

	double radius = 0;

	for (int i = 0; i < nFacets; ++i)
	{
		for (int j = 0; j < facets[i].size; ++j)
		{
			Point3f v_len = facets[i].arr[j] - p0;
			double len = Length(v_len);

			if (len > radius)
			{
				radius = len;
			}
		}
	}

	return radius;
}

double Particle::Area()
{
	double area = 0;

	for (int i = 0; i < nFacets; ++i)
	{
		area += facets[i].Area();
	}

	return area;
}

const complex &Particle::GetRefractiveIndex() const
{
	return m_refractiveIndex;
}

const Symmetry &Particle::GetSymmetry() const
{
	return m_symmetry;
}

void Particle::Move(float dx, float dy, float dz)
{
	for (int i = 0; i < nFacets; ++i)
	{
		for (int j = 0; j < defaultFacets[i].size; ++j)
		{
			facets[i].arr[j] = defaultFacets[i].arr[j] + Point3f(dx, dy, dz);
		}
	}
}

bool Particle::IsConcave() const
{
	return isConcave;
}

void Particle::Output()
{
	std::ofstream M("particle.dat", std::ios::out);

	for (int i = 0; i < nFacets; ++i)
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

void Particle::SetRefractiveIndex(const complex &value)
{
	m_refractiveIndex = value;
}

void Particle::SetDefaultNormals()
{
	for (int i = 0; i < nFacets; ++i)
	{
		defaultFacets[i].SetNormal();
	}
}

void Particle::Reset()
{
	for (int i = 0; i < nFacets; ++i)
	{
		facets[i] = defaultFacets[i];
	}
}

void Particle::SetDefaultCenters()
{
	for (int i = 0; i < nFacets; ++i)
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
	for (int i = 0; i < nFacets; ++i)
	{
		RotatePoint(defaultFacets[i].in_normal, facets[i].in_normal);
	}

	SetDParams();

	for (int i = 0; i < nFacets; ++i)
	{
		facets[i].ex_normal = -facets[i].in_normal;
		facets[i].ex_normal.d_param = -facets[i].in_normal.d_param;
	}
}

void Particle::SetDParams()
{
	for (int i = 0; i < nFacets; ++i)
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
