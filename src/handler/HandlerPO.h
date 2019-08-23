#pragma once

#include "Handler.h"

<<<<<<< HEAD
class HandlerPO : public Handler
{
public:
	HandlerPO(Particle *particle, Light *incidentLight, double wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
	void WriteMatricesToFile(std::string &destName) override;

	void SetScatteringSphere(const ScatteringSphere &grid);

	Arr2D M;				// Mueller matrices

protected:
	virtual void AddToMueller();

	void ComputeOpticalLengths(const Beam &beam, BeamInfo &info);

	void RotateJones(const Beam &beam, const BeamInfo &info,
					 const Vector3d &vf, const Vector3d &direction,
					 matrixC &matrix) const;
	void CleanJ();
	matrixC ComputeFnJones(const Matrix2x2c &matrix, const BeamInfo &info,
						   const Vector3d &direction);

	matrixC ApplyDiffraction(const Beam &beam, const BeamInfo &info,
							 const Vector3d &direction, const Vector3d &vf);

	BeamInfo ComputeBeamInfo(const Beam &beam);


protected:
	std::vector<Arr2DC> m_diffractedMatrices;	// Jones matrices
	bool isNanOccured = false;
	bool isNan = false;
};

=======

class HandlerPO : public Handler
{
public:
    HandlerPO(Particle *particle, Light *incidentLight, float wavelength = 0);

    void HandleBeams(std::vector<Beam> &beams) override;
    void WriteMatricesToFile(std::string &destName) override;

    void SetScatteringConus(const Conus &conus);

    void setCon20(bool value);

protected:
    void ApplyDiffraction(const Beam &beam, const Point3f &beamBasis,
                          const Vector3d &vf, const Vector3d &vr,
                          const matrixC &fnJones, matrixC &jones);

    void RotateJones(const Beam &beam, const Vector3f &T,
                     const Vector3d &vf, const Vector3d &vr, matrixC &J);

    void CleanJ();
    void AddToMueller();
    matrixC ComputeFnJones(const Matrix2x2c &jones, const Point3d &center,
                           const Vector3d &vr, double projLenght);

protected:
    std::vector<Arr2DC> J;	// Jones matrices
    Arr2D M;				// Mueller matrices

    Conus m_conus;			// back scattering conus
    bool isNanOccured = false;
    bool isNan = false;
    bool con20 = true;
};
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
