#ifndef __PETSC_PARALLEL_MANAGER_TURBULENT_KEPS_H__
#define __PETSC_PARALLEL_MANAGER_TURBULENT_KEPS_H__

#include "PetscParallelManagerTurbulent.h"
#include "../stencils/KineticEnergyBufferFillStencil.h"
#include "../stencils/KineticEnergyBufferReadStencil.h"
#include "../stencils/DissipationRateBufferFillStencil.h"
#include "../stencils/DissipationRateBufferReadStencil.h"

class PetscParallelManagerTurbulentKeps: public PetscParallelManagerTurbulent {

	private:
	TurbulentFlowField &_turbulentFlowField;

	FLOAT *_kineticEnergySendBufferLeftWall;
	FLOAT *_kineticEnergyRecvBufferLeftWall;
	FLOAT *_kineticEnergySendBufferRightWall;
	FLOAT *_kineticEnergyRecvBufferRightWall;
	FLOAT *_kineticEnergySendBufferTopWall;
	FLOAT *_kineticEnergyRecvBufferTopWall;
	FLOAT *_kineticEnergySendBufferBottomWall;
	FLOAT *_kineticEnergyRecvBufferBottomWall;
	FLOAT *_kineticEnergySendBufferFrontWall;
	FLOAT *_kineticEnergyRecvBufferFrontWall;
	FLOAT *_kineticEnergySendBufferBackWall;
	FLOAT *_kineticEnergyRecvBufferBackWall;
	
	FLOAT *_dissipationRateSendBufferLeftWall;
	FLOAT *_dissipationRateRecvBufferLeftWall;
	FLOAT *_dissipationRateSendBufferRightWall;
	FLOAT *_dissipationRateRecvBufferRightWall;
	FLOAT *_dissipationRateSendBufferTopWall;
	FLOAT *_dissipationRateRecvBufferTopWall;
	FLOAT *_dissipationRateSendBufferBottomWall;
	FLOAT *_dissipationRateRecvBufferBottomWall;
	FLOAT *_dissipationRateSendBufferFrontWall;
	FLOAT *_dissipationRateRecvBufferFrontWall;
	FLOAT *_dissipationRateSendBufferBackWall;
	FLOAT *_dissipationRateRecvBufferBackWall;

	KineticEnergyBufferFillStencil _kineticEnergyBufferFillStencil;
	KineticEnergyBufferReadStencil _kineticEnergyBufferReadStencil;
	DissipationRateBufferFillStencil _dissipationRateBufferFillStencil;
	DissipationRateBufferReadStencil _dissipationRateBufferReadStencil;

	ParallelBoundaryIterator<TurbulentFlowField> _kineticEnergyBufferFillIterator;
	ParallelBoundaryIterator<TurbulentFlowField> _kineticEnergyBufferReadIterator;
	ParallelBoundaryIterator<TurbulentFlowField> _dissipationRateBufferFillIterator;
	ParallelBoundaryIterator<TurbulentFlowField> _dissipationRateBufferReadIterator;

	public:

	PetscParallelManagerTurbulentKeps(TurbulentFlowField &flowField, Parameters &parameters);
	~PetscParallelManagerTurbulentKeps();
	void communicateKineticEnergy();
	void communicateDissipationRate();
  
};

#endif
