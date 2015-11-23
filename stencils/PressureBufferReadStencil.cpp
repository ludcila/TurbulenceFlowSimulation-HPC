#include "PressureBufferReadStencil.h"

PressureBufferReadStencil::PressureBufferReadStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall): BoundaryStencil(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
}

PressureBufferReadStencil::PressureBufferReadStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil(parameters) {
}

void PressureBufferReadStencil::applyLeftWall(FlowField & flowField, int i, int j) {
	flowField.getPressure().getScalar(i, j) = _bufferLeftWall[j];
}

void PressureBufferReadStencil::applyRightWall(FlowField & flowField, int i, int j) {
	flowField.getPressure().getScalar(i, j) = _bufferRightWall[j];
}

void PressureBufferReadStencil::applyBottomWall(FlowField & flowField, int i, int j) {

}

void PressureBufferReadStencil::applyTopWall(FlowField & flowField, int i, int j) {

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

