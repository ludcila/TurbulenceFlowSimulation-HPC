#include "PressureBufferReadStencil.h"

PressureBufferReadStencil::PressureBufferReadStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil<FlowField>(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

void PressureBufferReadStencil::applyLeftWall(FlowField & flowField, int i, int j) {
	flowField.getPressure().getScalar(i+1, j) = _bufferLeftWall[j];
}

void PressureBufferReadStencil::applyRightWall(FlowField & flowField, int i, int j) {
	flowField.getPressure().getScalar(i, j) = _bufferRightWall[j];
}

void PressureBufferReadStencil::applyBottomWall(FlowField & flowField, int i, int j) {
	flowField.getPressure().getScalar(i, j+1) = _bufferBottomWall[i];
}

void PressureBufferReadStencil::applyTopWall(FlowField & flowField, int i, int j) {
	flowField.getPressure().getScalar(i, j) = _bufferTopWall[i];
}

void PressureBufferReadStencil::applyLeftWall(FlowField & flowField, int i, int j, int k) {

}

void PressureBufferReadStencil::applyRightWall(FlowField & flowField, int i, int j, int k) {

}

void PressureBufferReadStencil::applyBottomWall(FlowField & flowField, int i, int j, int k) {

}

void PressureBufferReadStencil::applyTopWall(FlowField & flowField, int i, int j, int k) {

}

void PressureBufferReadStencil::applyFrontWall(FlowField & flowField, int i, int j, int k) {

}

void PressureBufferReadStencil::applyBackWall(FlowField & flowField, int i, int j, int k) {

}

