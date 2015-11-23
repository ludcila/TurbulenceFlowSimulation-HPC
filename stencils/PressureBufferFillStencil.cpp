#include "PressureBufferFillStencil.h"

PressureBufferFillStencil::PressureBufferFillStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

void PressureBufferFillStencil::applyLeftWall(FlowField & flowField, int i, int j) {
	this->_bufferLeftWall[j] = flowField.getPressure().getScalar(i+2, j);
}

void PressureBufferFillStencil::applyRightWall(FlowField & flowField, int i, int j) {
	this->_bufferRightWall[j] = flowField.getPressure().getScalar(i-1, j);
}

void PressureBufferFillStencil::applyBottomWall(FlowField & flowField, int i, int j) {
	this->_bufferBottomWall[i] = flowField.getPressure().getScalar(i, j+2);
}

void PressureBufferFillStencil::applyTopWall(FlowField & flowField, int i, int j) {
	this->_bufferTopWall[i] = flowField.getPressure().getScalar(i, j-1);
}

void PressureBufferFillStencil::applyLeftWall(FlowField & flowField, int i, int j, int k) {

}

void PressureBufferFillStencil::applyRightWall(FlowField & flowField, int i, int j, int k) {

}

void PressureBufferFillStencil::applyBottomWall(FlowField & flowField, int i, int j, int k) {

}

void PressureBufferFillStencil::applyTopWall(FlowField & flowField, int i, int j, int k) {

}

void PressureBufferFillStencil::applyFrontWall(FlowField & flowField, int i, int j, int k) {

}

void PressureBufferFillStencil::applyBackWall(FlowField & flowField, int i, int j, int k) {

}

