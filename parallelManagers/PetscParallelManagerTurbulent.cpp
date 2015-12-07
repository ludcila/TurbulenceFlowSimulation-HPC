#include "PetscParallelManagerTurbulent.h"

PetscParallelManagerTurbulent::PetscParallelManagerTurbulent(TurbulentFlowField &flowField, Parameters &parameters) :
  //parents' constructor
  PetscParallelManager(flowField, parameters),

  //Buffers
  _viscositySendBufferLeftWall		(new FLOAT[_cellsLeftRight]),
  _viscosityRecvBufferLeftWall		(new FLOAT[_cellsLeftRight]),
  _viscositySendBufferRightWall	  (new FLOAT[_cellsLeftRight]),
  _viscosityRecvBufferRightWall	  (new FLOAT[_cellsLeftRight]),
  _viscositySendBufferTopWall		  (new FLOAT[_cellsTopBottom]),
  _viscosityRecvBufferTopWall	  	(new FLOAT[_cellsTopBottom]),
  _viscositySendBufferBottomWall	(new FLOAT[_cellsTopBottom]),
  _viscosityRecvBufferBottomWall	(new FLOAT[_cellsTopBottom]),
  _viscositySendBufferFrontWall	  (parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
  _viscosityRecvBufferFrontWall	  (parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
  _viscositySendBufferBackWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
  _viscosityRecvBufferBackWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),

  //Stencils
  _viscosityBufferFillStencil(parameters, _viscositySendBufferLeftWall, _viscositySendBufferRightWall, _viscositySendBufferTopWall, _viscositySendBufferBottomWall, _viscositySendBufferFrontWall, _viscositySendBufferBackWall),
  _viscosityBufferReadStencil(parameters, _viscosityRecvBufferLeftWall, _viscosityRecvBufferRightWall, _viscosityRecvBufferTopWall, _viscosityRecvBufferBottomWall, _viscosityRecvBufferFrontWall, _viscosityRecvBufferBackWall),

  _viscosityBufferFillIterator(flowField, parameters, _viscosityBufferFillStencil),
  _viscosityBufferReadIterator(flowField, parameters, _viscosityBufferReadStencil)

{

}

PetscParallelManagerTurbulent::~PetscParallelManagerTurbulent(){
  delete [] _viscositySendBufferLeftWall;
  delete [] _viscosityRecvBufferLeftWall;
  delete [] _viscositySendBufferRightWall;
  delete [] _viscosityRecvBufferRightWall;
  delete [] _viscositySendBufferTopWall;
  delete [] _viscosityRecvBufferTopWall;
  delete [] _viscositySendBufferBottomWall;
  delete [] _viscositySendBufferFrontWall;
  delete [] _viscosityRecvBufferFrontWall;
  delete [] _viscositySendBufferBackWall;
  delete [] _viscosityRecvBufferBackWall;
}

void PetscParallelManagerTurbulent::communicateViscosity() {

	_viscosityBufferFillIterator.iterate();

	// Left to right & Right to left
	sendReceive(_viscositySendBufferRightWall, _parameters.parallel.rightNb, _viscosityRecvBufferLeftWall, _parameters.parallel.leftNb, _cellsLeftRight);
	sendReceive(_viscositySendBufferLeftWall, _parameters.parallel.leftNb, _viscosityRecvBufferRightWall, _parameters.parallel.rightNb, _cellsLeftRight);
	// Top to bottom & Bottom to top
	sendReceive(_viscositySendBufferTopWall, _parameters.parallel.topNb, _viscosityRecvBufferBottomWall, _parameters.parallel.bottomNb, _cellsTopBottom);
	sendReceive(_viscositySendBufferBottomWall, _parameters.parallel.bottomNb, _viscosityRecvBufferTopWall, _parameters.parallel.topNb, _cellsTopBottom);
	// Front to back & Back to front
	sendReceive(_viscositySendBufferFrontWall, _parameters.parallel.frontNb, _viscosityRecvBufferBackWall, _parameters.parallel.backNb, _cellsFrontBack);
	sendReceive(_viscositySendBufferBackWall, _parameters.parallel.backNb, _viscosityRecvBufferFrontWall, _parameters.parallel.frontNb, _cellsFrontBack);

	_viscosityBufferReadIterator.iterate();

}
