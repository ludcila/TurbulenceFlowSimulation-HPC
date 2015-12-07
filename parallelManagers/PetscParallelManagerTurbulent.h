#ifndef __PETSC_PARALLEL_MANAGER_TURBULENT_H__
#define __PETSC_PARALLEL_MANAGER_TURBULENT_H__

#include "PetscParallelManager.h"
#include "../stencils/ViscosityBufferFillStencil.h"
#include "../stencils/ViscosityBufferReadStencil.h"


class PetscParallelManagerTurbulent: public PetscParallelManager {

  FLOAT *_viscositySendBufferLeftWall;
  FLOAT *_viscosityRecvBufferLeftWall;
  FLOAT *_viscositySendBufferRightWall;
  FLOAT *_viscosityRecvBufferRightWall;
  FLOAT *_viscositySendBufferTopWall;
  FLOAT *_viscosityRecvBufferTopWall;
  FLOAT *_viscositySendBufferBottomWall;
  FLOAT *_viscosityRecvBufferBottomWall;
  FLOAT *_viscositySendBufferFrontWall;
  FLOAT *_viscosityRecvBufferFrontWall;
  FLOAT *_viscositySendBufferBackWall;
  FLOAT *_viscosityRecvBufferBackWall;

  ViscosityBufferFillStencil _viscosityBufferFillStencil;
	ViscosityBufferReadStencil _viscosityBufferReadStencil;

  ParallelBoundaryIterator<TurbulentFlowField> _viscosityBufferFillIterator;
  ParallelBoundaryIterator<TurbulentFlowField> _viscosityBufferReadIterator;

  PetscParallelManagerTurbulent(TurbulentFlowField &flowField, Parameters &parameters);
  ~PetscParallelManagerTurbulent();
  void communicateViscosity();

};

#endif
