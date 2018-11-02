#include "Rotator.h"

void LocalRotator::SetRotationAngle(const Angle3d &angle)
{
	double cosA, cosB, cosG,
			sinA, sinB, sinG;

	sincos(angle.alpha, &sinA, &cosA);
	sincos(angle.beta,  &sinB, &cosB);
	sincos(angle.gamma, &sinG, &cosG);

	double cosAcosB = cosA*cosB;
	double sinAcosG = sinA*cosG;
	double sinAsinG = sinA*sinG;

	m_rotationMatrix[0][0] = cosAcosB*cosG - sinAsinG;
	m_rotationMatrix[1][0] = sinAcosG*cosB + cosA*sinG;
	m_rotationMatrix[2][0] = -sinB*cosG;

	m_rotationMatrix[0][1] = -(cosAcosB*sinG + sinAcosG);
	m_rotationMatrix[1][1] = cosA*cosG - sinAsinG*cosB;
	m_rotationMatrix[2][1] = sinB*sinG;

	m_rotationMatrix[0][2] = cosA*sinB;
	m_rotationMatrix[1][2] = sinA*sinB;
	m_rotationMatrix[2][2] = cosB;
}

void GlobalRotator::SetRotationAngle(const Angle3d &angle)
{
	double cosF, cosT, cosP,
			sinF, sinT, sinP;

	sincos(angle.alpha, &sinP, &cosP);
	sincos(angle.beta, &sinF, &cosF);
	sincos(angle.gamma, &sinT, &cosT);

	double sinFsinT = sinF*sinT;
	double cosFsinT = cosF*sinT;

	m_rotationMatrix[0][0] = cosT*cosP;
	m_rotationMatrix[1][0] = sinFsinT*cosP + cosF*sinP;
	m_rotationMatrix[2][0] = -cosFsinT*cosP + sinF*sinP;

	m_rotationMatrix[0][1] = -cosT*sinP;
	m_rotationMatrix[1][1] = -sinFsinT*sinP + cosF*cosP;
	m_rotationMatrix[2][1] = cosFsinT*sinP + sinF*cosP;

	m_rotationMatrix[0][2] = sinT;
	m_rotationMatrix[1][2] = -sinF*cosT;
	m_rotationMatrix[2][2] = cosF*cosT;
}

void Rotator::RotatePoint(const Point3f &point, Point3f &result)
{
	result.coordinates[0] = point.coordinates[0] * m_rotationMatrix[0][0] +
			point.coordinates[1] * m_rotationMatrix[0][1] +
			point.coordinates[2] * m_rotationMatrix[0][2];

	result.coordinates[1] = point.coordinates[0] * m_rotationMatrix[1][0] +
			point.coordinates[1] * m_rotationMatrix[1][1] +
			point.coordinates[2] * m_rotationMatrix[1][2];

	result.coordinates[2] = point.coordinates[0] * m_rotationMatrix[2][0] +
			point.coordinates[1] * m_rotationMatrix[2][1] +
			point.coordinates[2] * m_rotationMatrix[2][2];
}
