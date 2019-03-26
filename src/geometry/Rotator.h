#pragma once

#include <geometry_lib.h>

/**
 * @brief Rotates points by given angle
 */
class Rotator
{
public:
    virtual void SetRotationAngle(const Orientation &angle) = 0;
    void RotatePoint(const Point3f &point, Point3f &result);

protected:
    double m_rotationMatrix[ROT_MTR_RANK][ROT_MTR_RANK];
};

class LocalRotator : public Rotator
{
public:
    void SetRotationAngle(const Orientation &angle);
};

class GlobalRotator : public Rotator
{
public:
    void SetRotationAngle(const Orientation &angle) override;
};
