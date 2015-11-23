#include "PetscParallelManager.h"
#include "../stencils/PressureBufferFillStencil.h"
#include "../stencils/PressureBufferReadStencil.h"

PetscParallelManager::PetscParallelManager(FlowField &flowField, Parameters &parameters) :
	
	_parameters(parameters),
	_flowField(flowField),
	
	// Dimensions
	_cellsX(flowField.getCellsX()),
	_cellsY(flowField.getCellsY()),
	_cellsZ(parameters.geometry.dim == 2 ? 1 : flowField.getCellsZ()),
	_cellsFrontBack(_cellsX * _cellsY),
	_cellsTopBottom(_cellsX * _cellsZ),
	_cellsLeftRight(_cellsY * _cellsZ),
	
	// Buffers
	_pressureSendBufferLeftWall		(new FLOAT[_cellsLeftRight]),
	_pressureRecvBufferLeftWall		(new FLOAT[_cellsLeftRight]),
	_pressureSendBufferRightWall	(new FLOAT[_cellsLeftRight]),
	_pressureRecvBufferRightWall	(new FLOAT[_cellsLeftRight]),
	_pressureSendBufferTopWall		(new FLOAT[_cellsTopBottom]),
	_pressureRecvBufferTopWall		(new FLOAT[_cellsTopBottom]),
	_pressureSendBufferBottomWall	(new FLOAT[_cellsTopBottom]),
	_pressureRecvBufferBottomWall	(new FLOAT[_cellsTopBottom]),
	_pressureSendBufferFrontWall	(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
	_pressureRecvBufferFrontWall	(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
	_pressureSendBufferBackWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
	_pressureRecvBufferBackWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
	_velocitySendBufferLeftWall		(new FLOAT[_cellsLeftRight * parameters.geometry.dim]),
	_velocityRecvBufferLeftWall		(new FLOAT[_cellsLeftRight * parameters.geometry.dim]),
	_velocitySendBufferRightWall	(new FLOAT[_cellsLeftRight * parameters.geometry.dim]),
	_velocityRecvBufferRightWall	(new FLOAT[_cellsLeftRight * parameters.geometry.dim]),
	_velocitySendBufferTopWall		(new FLOAT[_cellsTopBottom * parameters.geometry.dim]),
	_velocityRecvBufferTopWall		(new FLOAT[_cellsTopBottom * parameters.geometry.dim]),
	_velocitySendBufferBottomWall	(new FLOAT[_cellsTopBottom * parameters.geometry.dim]),
	_velocityRecvBufferBottomWall	(new FLOAT[_cellsTopBottom * parameters.geometry.dim]),
	_velocitySendBufferFrontWall	(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack * parameters.geometry.dim]),
	_velocityRecvBufferFrontWall	(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack * parameters.geometry.dim]),
	_velocitySendBufferBackWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack * parameters.geometry.dim]),
	_velocityRecvBufferBackWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack * parameters.geometry.dim]),
	
	// Stencils
	_pressureBufferFillStencil(parameters, _pressureSendBufferLeftWall, _pressureSendBufferRightWall, _pressureSendBufferTopWall, _pressureSendBufferBottomWall, _pressureSendBufferFrontWall, _pressureSendBufferBackWall),
	_pressureBufferReadStencil(parameters, _pressureRecvBufferLeftWall, _pressureRecvBufferRightWall, _pressureRecvBufferTopWall, _pressureRecvBufferBottomWall, _pressureRecvBufferFrontWall, _pressureRecvBufferBackWall),
	// _velocityBufferFillStencil(parameters, _velocitySendBufferLeftWall, _velocitySendBufferRightWall, _velocitySendBufferTopWall, _velocitySendBufferBottomWall, _velocitySendBufferFrontWall, _velocitySendBufferBackWall),
	// _velocityBufferReadStencil(parameters, _velocityRecvBufferLeftWall, _velocityRecvBufferRightWall, _velocityRecvBufferTopWall, _velocityRecvBufferBottomWall, _velocityRecvBufferFrontWall, _velocityRecvBufferBackWall),

	// Iterators
	_pressureBufferFillIterator(flowField, parameters, _pressureBufferFillStencil),
	_pressureBufferReadIterator(flowField, parameters, _pressureBufferReadStencil)
	// _velocityBufferFillIterator(flowField, parameters, _velocityBufferFillStencil),
	// _velocityBufferReadIterator(flowField, parameters, _velocityBufferReadStencil)
	
{

}

PetscParallelManager::~PetscParallelManager() {
	delete [] _pressureSendBufferLeftWall;
	delete [] _pressureRecvBufferLeftWall;
	delete [] _pressureSendBufferRightWall;
	delete [] _pressureRecvBufferRightWall;
	delete [] _pressureSendBufferTopWall;
	delete [] _pressureRecvBufferTopWall;
	delete [] _pressureSendBufferBottomWall;
	delete [] _pressureSendBufferFrontWall;
	delete [] _pressureRecvBufferFrontWall;
	delete [] _pressureSendBufferBackWall;
	delete [] _pressureRecvBufferBackWall;
	delete [] _velocitySendBufferLeftWall;
	delete [] _velocityRecvBufferLeftWall;
	delete [] _velocitySendBufferRightWall;
	delete [] _velocityRecvBufferRightWall;
	delete [] _velocitySendBufferTopWall;
	delete [] _velocityRecvBufferTopWall;
	delete [] _velocitySendBufferBottomWall;
	delete [] _velocityRecvBufferBottomWall;
	delete [] _velocitySendBufferFrontWall;
	delete [] _velocityRecvBufferFrontWall;
	delete [] _velocitySendBufferBackWall;
	delete [] _velocityRecvBufferBackWall;
}

void PetscParallelManager::communicatePressure() {
	
	_pressureBufferFillIterator.iterate();
	
	// Left to right & Right to left
	sendReceive(_pressureSendBufferRightWall, _parameters.parallel.rightNb, _pressureRecvBufferLeftWall, _parameters.parallel.leftNb, _cellsLeftRight);
	sendReceive(_pressureSendBufferLeftWall, _parameters.parallel.leftNb, _pressureRecvBufferRightWall, _parameters.parallel.rightNb, _cellsLeftRight);
	
	_pressureBufferReadIterator.iterate();
	
}

void PetscParallelManager::sendReceive(FLOAT *sendBuffer, int sendTo, FLOAT *receiveBuffer, int receiveFrom, int size) {
	MPI_Sendrecv(
		sendBuffer, size, MPI_DOUBLE, sendTo, 0,
		receiveBuffer, size, MPI_DOUBLE, receiveFrom, 0,
		MPI_COMM_WORLD, MPI_STATUS_IGNORE
	);
	// MPI_DOUBLE ok?
}

