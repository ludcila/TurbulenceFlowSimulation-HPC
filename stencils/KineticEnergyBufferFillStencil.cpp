#include "KineticEnergyBufferFillStencil.h"

KineticEnergyBufferFillStencil::KineticEnergyBufferFillStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil<TurbulentFlowField>(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

// 2D stencils
void KineticEnergyBufferFillStencil::applyLeftWall(TurbulentFlowField & TurbulentFlowField, int i, int j) {
	this->_bufferLeftWall[j] = TurbulentFlowField.getKineticEnergy().getScalar(i+2, j);
}

void KineticEnergyBufferFillStencil::applyRightWall(TurbulentFlowField & TurbulentFlowField, int i, int j) {
	this->_bufferRightWall[j] = TurbulentFlowField.getKineticEnergy().getScalar(i-1, j);
}

void KineticEnergyBufferFillStencil::applyBottomWall(TurbulentFlowField & TurbulentFlowField, int i, int j) {
	this->_bufferBottomWall[i] = TurbulentFlowField.getKineticEnergy().getScalar(i, j+2);
}

void KineticEnergyBufferFillStencil::applyTopWall(TurbulentFlowField & TurbulentFlowField, int i, int j) {
	this->_bufferTopWall[i] = TurbulentFlowField.getKineticEnergy().getScalar(i, j-1);
}


// 3D stencils
void KineticEnergyBufferFillStencil::applyLeftWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	this->_bufferLeftWall[TurbulentFlowField.getCellsZ()*j+k] = TurbulentFlowField.getKineticEnergy().getScalar(i+2, j, k);
}

void KineticEnergyBufferFillStencil::applyRightWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	this->_bufferRightWall[TurbulentFlowField.getCellsZ()*j+k] = TurbulentFlowField.getKineticEnergy().getScalar(i-1, j, k);
}

void KineticEnergyBufferFillStencil::applyBottomWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	this->_bufferBottomWall[TurbulentFlowField.getCellsZ()*i+k] = TurbulentFlowField.getKineticEnergy().getScalar(i, j+2, k);
}

void KineticEnergyBufferFillStencil::applyTopWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	this->_bufferTopWall[TurbulentFlowField.getCellsZ()*i+k] = TurbulentFlowField.getKineticEnergy().getScalar(i, j-1, k);
}

void KineticEnergyBufferFillStencil::applyFrontWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	this->_bufferFrontWall[TurbulentFlowField.getCellsY()*i+j] = TurbulentFlowField.getKineticEnergy().getScalar(i, j, k+2);
}

void KineticEnergyBufferFillStencil::applyBackWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	this->_bufferBackWall[TurbulentFlowField.getCellsY()*i+j] = TurbulentFlowField.getKineticEnergy().getScalar(i, j, k-1);
}

