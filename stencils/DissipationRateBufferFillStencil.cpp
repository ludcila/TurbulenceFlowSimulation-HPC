#include "DissipationRateBufferFillStencil.h"

DissipationRateBufferFillStencil::DissipationRateBufferFillStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil<TurbulentFlowField>(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

// 2D stencils
void DissipationRateBufferFillStencil::applyLeftWall(TurbulentFlowField & TurbulentFlowField, int i, int j) {
	this->_bufferLeftWall[j] = TurbulentFlowField.getDissipationRate().getScalar(i+2, j);
}

void DissipationRateBufferFillStencil::applyRightWall(TurbulentFlowField & TurbulentFlowField, int i, int j) {
	this->_bufferRightWall[j] = TurbulentFlowField.getDissipationRate().getScalar(i-1, j);
}

void DissipationRateBufferFillStencil::applyBottomWall(TurbulentFlowField & TurbulentFlowField, int i, int j) {
	this->_bufferBottomWall[i] = TurbulentFlowField.getDissipationRate().getScalar(i, j+2);
}

void DissipationRateBufferFillStencil::applyTopWall(TurbulentFlowField & TurbulentFlowField, int i, int j) {
	this->_bufferTopWall[i] = TurbulentFlowField.getDissipationRate().getScalar(i, j-1);
}


// 3D stencils
void DissipationRateBufferFillStencil::applyLeftWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	this->_bufferLeftWall[TurbulentFlowField.getCellsZ()*j+k] = TurbulentFlowField.getDissipationRate().getScalar(i+2, j, k);
}

void DissipationRateBufferFillStencil::applyRightWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	this->_bufferRightWall[TurbulentFlowField.getCellsZ()*j+k] = TurbulentFlowField.getDissipationRate().getScalar(i-1, j, k);
}

void DissipationRateBufferFillStencil::applyBottomWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	this->_bufferBottomWall[TurbulentFlowField.getCellsZ()*i+k] = TurbulentFlowField.getDissipationRate().getScalar(i, j+2, k);
}

void DissipationRateBufferFillStencil::applyTopWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	this->_bufferTopWall[TurbulentFlowField.getCellsZ()*i+k] = TurbulentFlowField.getDissipationRate().getScalar(i, j-1, k);
}

void DissipationRateBufferFillStencil::applyFrontWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	this->_bufferFrontWall[TurbulentFlowField.getCellsY()*i+j] = TurbulentFlowField.getDissipationRate().getScalar(i, j, k+2);
}

void DissipationRateBufferFillStencil::applyBackWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k) {
	this->_bufferBackWall[TurbulentFlowField.getCellsY()*i+j] = TurbulentFlowField.getDissipationRate().getScalar(i, j, k-1);
}

