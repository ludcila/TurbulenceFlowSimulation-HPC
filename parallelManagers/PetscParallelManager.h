#ifndef __PETSC_PARALLEL_MANAGER_H__
#define __PETSC_PARALLEL_MANAGER_H__

#include "../FlowField.h"
#include "../Parameters.h"
#include "../Iterators.h"
#include "../stencils/PressureBufferFillStencil.h"
#include "../stencils/PressureBufferReadStencil.h"

class PetscParallelManager {

	protected:
	
	int _cellsX;
	int _cellsY;
	int _cellsZ; // Set to 1 for 2D case
	
	// Number of elements in the parallel boundaries
	int _cellsLeftRight;
	int _cellsTopBottom;
	int _cellsFrontBack;
	
	Parameters &_parameters;
	FlowField &_flowField;
	
	// Buffers
	FLOAT *_pressureSendBufferLeftWall = NULL;
	FLOAT *_pressureRecvBufferLeftWall = NULL;
	FLOAT *_pressureSendBufferRightWall = NULL;
	FLOAT *_pressureRecvBufferRightWall = NULL;
	FLOAT *_pressureSendBufferTopWall = NULL;
	FLOAT *_pressureRecvBufferTopWall = NULL;
	FLOAT *_pressureSendBufferBottomWall = NULL;
	FLOAT *_pressureRecvBufferBottomWall = NULL;
	FLOAT *_pressureSendBufferFrontWall = NULL;
	FLOAT *_pressureRecvBufferFrontWall = NULL;
	FLOAT *_pressureSendBufferBackWall = NULL;
	FLOAT *_pressureRecvBufferBackWall = NULL;
	FLOAT *_velocitySendBufferLeftWall = NULL;
	FLOAT *_velocityRecvBufferLeftWall = NULL;
	FLOAT *_velocitySendBufferRightWall = NULL;
	FLOAT *_velocityRecvBufferRightWall = NULL;
	FLOAT *_velocitySendBufferTopWall = NULL;
	FLOAT *_velocityRecvBufferTopWall = NULL;
	FLOAT *_velocitySendBufferBottomWall = NULL;
	FLOAT *_velocityRecvBufferBottomWall = NULL;
	FLOAT *_velocitySendBufferFrontWall = NULL;
	FLOAT *_velocityRecvBufferFrontWall = NULL;
	FLOAT *_velocitySendBufferBackWall = NULL;
	FLOAT *_velocityRecvBufferBackWall = NULL;
	
	// Stencils
	PressureBufferFillStencil _pressureBufferFillStencil;
	PressureBufferReadStencil _pressureBufferReadStencil;
	// VelocityBufferFillStencil _velocityBufferFillStencil;
	// VelocityBufferReadStencil _velocityBufferReadStencil;
	
	// Iterators
	ParallelBoundaryIterator<FlowField> _pressureBufferFillIterator;
	ParallelBoundaryIterator<FlowField> _pressureBufferReadIterator;
	//ParallelBoundaryIterator _velocityBufferFillIterator;
	//ParallelBoundaryIterator _velocityBufferReadIterator;

	public:
	
	PetscParallelManager(FlowField &flowField, Parameters &parameters);
	~PetscParallelManager();
	void communicatePressure();
	void sendReceive(FLOAT *sendBuffer, int sendTo, FLOAT *receiveBuffer, int receiveFrom, int size);
	

};

#endif
