#include "VelocityBufferReadStencil.h"

VelocityBufferReadStencil::VelocityBufferReadStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

void VelocityBufferReadStencil::applyLeftWall(FlowField & flowField, int i, int j) {
	flowField.getVelocity().getVector(i, j, 0)[0] = _bufferLeftWall[BUFFER_2D_POS_X(j)];
	flowField.getVelocity().getVector(i, j, 0)[1] = _bufferLeftWall[BUFFER_2D_POS_Y(j)];
}

void VelocityBufferReadStencil::applyRightWall(FlowField & flowField, int i, int j) {
	flowField.getVelocity().getVector(i, j, 0)[0] = _bufferRightWall[BUFFER_2D_POS_X(j)];
	flowField.getVelocity().getVector(i, j, 0)[1] = _bufferRightWall[BUFFER_2D_POS_Y(j)];
}

void VelocityBufferReadStencil::applyBottomWall(FlowField & flowField, int i, int j) {
	flowField.getVelocity().getVector(i, j, 0)[0] = _bufferBottomWall[BUFFER_2D_POS_X(i)];
	flowField.getVelocity().getVector(i, j, 0)[1] = _bufferBottomWall[BUFFER_2D_POS_Y(i)];
}

void VelocityBufferReadStencil::applyTopWall(FlowField & flowField, int i, int j) {
	flowField.getVelocity().getVector(i, j, 0)[0] = _bufferTopWall[BUFFER_2D_POS_X(i)];
	flowField.getVelocity().getVector(i, j, 0)[1] = _bufferTopWall[BUFFER_2D_POS_Y(i)];
}

void VelocityBufferReadStencil::applyLeftWall(FlowField & flowField, int i, int j, int k) {

}

void VelocityBufferReadStencil::applyRightWall(FlowField & flowField, int i, int j, int k) {

}

void VelocityBufferReadStencil::applyBottomWall(FlowField & flowField, int i, int j, int k) {

}

void VelocityBufferReadStencil::applyTopWall(FlowField & flowField, int i, int j, int k) {

}

void VelocityBufferReadStencil::applyFrontWall(FlowField & flowField, int i, int j, int k) {

}

void VelocityBufferReadStencil::applyBackWall(FlowField & flowField, int i, int j, int k) {

}

