#pragma once

#include "Handler.h"

class HandlerPO : public Handler
{
public:
	HandlerPO(Particle *particle, Light *incidentLight, double wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
	void WriteMatricesToFile(std::string &destName) override;

	void SetScatteringSphere(const ScatteringSphere &grid);

	Arr2D M;				// Mueller matrices

	matrixC ApplyDiffraction(const Beam &beam, const BeamInfo &info,
							 const Vector3d &direction, const Vector3d &vf);

protected:
	virtual void AddToMueller();

	void RotateJones(const Beam &beam, const BeamInfo &info,
					 const Vector3d &vf, const Vector3d &direction,
					 matrixC &matrix) const;
	void CleanJ();
	matrixC ComputeFnJones(const Matrix2x2c &matrix, const BeamInfo &info,
						   const Vector3d &direction);


protected:
	std::vector<Arr2DC> m_diffractedMatrices;	// Jones matrices
	bool isNanOccured = false;
	bool isNan = false;
};
