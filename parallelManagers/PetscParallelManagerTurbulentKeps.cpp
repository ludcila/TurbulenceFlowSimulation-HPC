#include "PetscParallelManagerTurbulentKeps.h"
#include "PetscParallelConfiguration.h"

PetscParallelManagerTurbulentKeps::PetscParallelManagerTurbulentKeps(TurbulentFlowField &flowField, Parameters &parameters) :

	//parents' constructor
	PetscParallelManagerTurbulent(flowField, parameters),
	_turbulentFlowField(flowField),

	//Buffers
	_kineticEnergySendBufferLeftWall		(new FLOAT[_cellsLeftRight]),
	_kineticEnergyRecvBufferLeftWall		(new FLOAT[_cellsLeftRight]),
	_kineticEnergySendBufferRightWall	  	(new FLOAT[_cellsLeftRight]),
	_kineticEnergyRecvBufferRightWall	  	(new FLOAT[_cellsLeftRight]),
	_kineticEnergySendBufferTopWall		  	(new FLOAT[_cellsTopBottom]),
	_kineticEnergyRecvBufferTopWall	  		(new FLOAT[_cellsTopBottom]),
	_kineticEnergySendBufferBottomWall		(new FLOAT[_cellsTopBottom]),
	_kineticEnergyRecvBufferBottomWall		(new FLOAT[_cellsTopBottom]),
	_kineticEnergySendBufferFrontWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
	_kineticEnergyRecvBufferFrontWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
	_kineticEnergySendBufferBackWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
	_kineticEnergyRecvBufferBackWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
	
	_dissipationRateSendBufferLeftWall		(new FLOAT[_cellsLeftRight]),
	_dissipationRateRecvBufferLeftWall		(new FLOAT[_cellsLeftRight]),
	_dissipationRateSendBufferRightWall	  	(new FLOAT[_cellsLeftRight]),
	_dissipationRateRecvBufferRightWall	  	(new FLOAT[_cellsLeftRight]),
	_dissipationRateSendBufferTopWall		(new FLOAT[_cellsTopBottom]),
	_dissipationRateRecvBufferTopWall	  	(new FLOAT[_cellsTopBottom]),
	_dissipationRateSendBufferBottomWall	(new FLOAT[_cellsTopBottom]),
	_dissipationRateRecvBufferBottomWall	(new FLOAT[_cellsTopBottom]),
	_dissipationRateSendBufferFrontWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
	_dissipationRateRecvBufferFrontWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
	_dissipationRateSendBufferBackWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
	_dissipationRateRecvBufferBackWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),

	//Stencils
	_kineticEnergyBufferFillStencil(parameters, _kineticEnergySendBufferLeftWall, _kineticEnergySendBufferRightWall, _kineticEnergySendBufferTopWall, _kineticEnergySendBufferBottomWall, _kineticEnergySendBufferFrontWall, _kineticEnergySendBufferBackWall),
	_kineticEnergyBufferReadStencil(parameters, _kineticEnergyRecvBufferLeftWall, _kineticEnergyRecvBufferRightWall, _kineticEnergyRecvBufferTopWall, _kineticEnergyRecvBufferBottomWall, _kineticEnergyRecvBufferFrontWall, _kineticEnergyRecvBufferBackWall),
	
	_dissipationRateBufferFillStencil(parameters, _dissipationRateSendBufferLeftWall, _dissipationRateSendBufferRightWall, _dissipationRateSendBufferTopWall, _dissipationRateSendBufferBottomWall, _dissipationRateSendBufferFrontWall, _dissipationRateSendBufferBackWall),
	_dissipationRateBufferReadStencil(parameters, _dissipationRateRecvBufferLeftWall, _dissipationRateRecvBufferRightWall, _dissipationRateRecvBufferTopWall, _dissipationRateRecvBufferBottomWall, _dissipationRateRecvBufferFrontWall, _dissipationRateRecvBufferBackWall),

	//Iterators
	_kineticEnergyBufferFillIterator(flowField, parameters, _kineticEnergyBufferFillStencil),
	_kineticEnergyBufferReadIterator(flowField, parameters, _kineticEnergyBufferReadStencil),
	
	_dissipationRateBufferFillIterator(flowField, parameters, _dissipationRateBufferFillStencil),
	_dissipationRateBufferReadIterator(flowField, parameters, _dissipationRateBufferReadStencil)

{

}

PetscParallelManagerTurbulentKeps::~PetscParallelManagerTurbulentKeps(){
	delete [] _kineticEnergySendBufferLeftWall;
	delete [] _kineticEnergyRecvBufferLeftWall;
	delete [] _kineticEnergySendBufferRightWall;
	delete [] _kineticEnergyRecvBufferRightWall;
	delete [] _kineticEnergySendBufferTopWall;
	delete [] _kineticEnergyRecvBufferTopWall;
	delete [] _kineticEnergySendBufferBottomWall;
	delete [] _kineticEnergyRecvBufferBottomWall;
	delete [] _kineticEnergySendBufferFrontWall;
	delete [] _kineticEnergyRecvBufferFrontWall;
	delete [] _kineticEnergySendBufferBackWall;
	delete [] _kineticEnergyRecvBufferBackWall;
	delete [] _dissipationRateSendBufferLeftWall;
	delete [] _dissipationRateRecvBufferLeftWall;
	delete [] _dissipationRateSendBufferRightWall;
	delete [] _dissipationRateRecvBufferRightWall;
	delete [] _dissipationRateSendBufferTopWall;
	delete [] _dissipationRateRecvBufferTopWall;
	delete [] _dissipationRateSendBufferBottomWall;
	delete [] _dissipationRateRecvBufferBottomWall;
	delete [] _dissipationRateSendBufferFrontWall;
	delete [] _dissipationRateRecvBufferFrontWall;
	delete [] _dissipationRateSendBufferBackWall;
	delete [] _dissipationRateRecvBufferBackWall;
}

void PetscParallelManagerTurbulentKeps::communicateKineticEnergy() {

	_kineticEnergyBufferFillIterator.iterate();

	// Left to right & Right to left
	sendReceive(_kineticEnergySendBufferRightWall, _parameters.parallel.rightNb, _kineticEnergyRecvBufferLeftWall, _parameters.parallel.leftNb, _cellsLeftRight);
	sendReceive(_kineticEnergySendBufferLeftWall, _parameters.parallel.leftNb, _kineticEnergyRecvBufferRightWall, _parameters.parallel.rightNb, _cellsLeftRight);
	// Top to bottom & Bottom to top
	sendReceive(_kineticEnergySendBufferTopWall, _parameters.parallel.topNb, _kineticEnergyRecvBufferBottomWall, _parameters.parallel.bottomNb, _cellsTopBottom);
	sendReceive(_kineticEnergySendBufferBottomWall, _parameters.parallel.bottomNb, _kineticEnergyRecvBufferTopWall, _parameters.parallel.topNb, _cellsTopBottom);
	// Front to back & Back to front
	sendReceive(_kineticEnergySendBufferFrontWall, _parameters.parallel.frontNb, _kineticEnergyRecvBufferBackWall, _parameters.parallel.backNb, _cellsFrontBack);
	sendReceive(_kineticEnergySendBufferBackWall, _parameters.parallel.backNb, _kineticEnergyRecvBufferFrontWall, _parameters.parallel.frontNb, _cellsFrontBack);

	_kineticEnergyBufferReadIterator.iterate();

}

void PetscParallelManagerTurbulentKeps::communicateDissipationRate() {

	_dissipationRateBufferFillIterator.iterate();

	// Left to right & Right to left
	sendReceive(_dissipationRateSendBufferRightWall, _parameters.parallel.rightNb, _dissipationRateRecvBufferLeftWall, _parameters.parallel.leftNb, _cellsLeftRight);
	sendReceive(_dissipationRateSendBufferLeftWall, _parameters.parallel.leftNb, _dissipationRateRecvBufferRightWall, _parameters.parallel.rightNb, _cellsLeftRight);
	// Top to bottom & Bottom to top
	sendReceive(_dissipationRateSendBufferTopWall, _parameters.parallel.topNb, _dissipationRateRecvBufferBottomWall, _parameters.parallel.bottomNb, _cellsTopBottom);
	sendReceive(_dissipationRateSendBufferBottomWall, _parameters.parallel.bottomNb, _dissipationRateRecvBufferTopWall, _parameters.parallel.topNb, _cellsTopBottom);
	// Front to back & Back to front
	sendReceive(_dissipationRateSendBufferFrontWall, _parameters.parallel.frontNb, _dissipationRateRecvBufferBackWall, _parameters.parallel.backNb, _cellsFrontBack);
	sendReceive(_dissipationRateSendBufferBackWall, _parameters.parallel.backNb, _dissipationRateRecvBufferFrontWall, _parameters.parallel.frontNb, _cellsFrontBack);

	_dissipationRateBufferReadIterator.iterate();

}


