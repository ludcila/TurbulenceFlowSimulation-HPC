#include "KineticEnergyBufferReadStencil.h"

KineticEnergyBufferReadStencil::KineticEnergyBufferReadStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil<TurbulentFlowField>(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

// 2D stencils
void KineticEnergyBufferReadStencil::applyLeftWall(TurbulentFlowField & TurbulentFlowField, int i, int j) {
	TurbulentFlowField.getKineticEnergy().getScalar(i+1, j) = _bufferLeftWall[j];
}

void KineticEnergyBufferReadStencil::applyRightWall(TurbulentFlowField & TurbulentFlowField, int i, int j) {
	TurbulentFlowField.getKineticEnergy().getScalar(i, j) = _bufferRightWall[j];
}

void KineticEnergyBufferReadStencil::applyBottomWall(TurbulentFlowField & TurbulentFlowField, int i, int j) {
	TurbulentFlowField.getKineticEnergy().getScalar(i, j+1) = _bufferBottomWall[i];
}

void KineticEnergyBufferReadStencil::applyTopWall(TurbulentFlowField & TurbulentFlowField, int i, int j) {
	TurbulentFlowField.getKineticEnergy().getScalar(i, j) = _bufferTopWall[i];
}



// 3D stencils
void KineticEnergyBufferReadStencil::applyLeftWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	TurbulentFlowField.getKineticEnergy().getScalar(i+1, j, k) = _bufferLeftWall[TurbulentFlowField.getCellsZ()*j+k];
}

void KineticEnergyBufferReadStencil::applyRightWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	TurbulentFlowField.getKineticEnergy().getScalar(i, j, k) = _bufferRightWall[TurbulentFlowField.getCellsZ()*j+k];
}

void KineticEnergyBufferReadStencil::applyBottomWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	TurbulentFlowField.getKineticEnergy().getScalar(i, j+1, k) = _bufferBottomWall[TurbulentFlowField.getCellsZ()*i+k];
}

void KineticEnergyBufferReadStencil::applyTopWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	TurbulentFlowField.getKineticEnergy().getScalar(i, j, k) = _bufferTopWall[TurbulentFlowField.getCellsZ()*i+k];
}

void KineticEnergyBufferReadStencil::applyFrontWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	TurbulentFlowField.getKineticEnergy().getScalar(i, j, k+1) = _bufferFrontWall[TurbulentFlowField.getCellsY()*i+j];
}

void KineticEnergyBufferReadStencil::applyBackWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	TurbulentFlowField.getKineticEnergy().getScalar(i, j, k) = _bufferBackWall[TurbulentFlowField.getCellsY()*i+j];
}

