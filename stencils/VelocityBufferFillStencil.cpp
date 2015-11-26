#include "VelocityBufferFillStencil.h"

VelocityBufferFillStencil::VelocityBufferFillStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

void VelocityBufferFillStencil::applyLeftWall(FlowField & flowField, int i, int j) {
	this->_bufferLeftWall[BUFFER_2D_POS_X(j)] = flowField.getVelocity().getVector(i+2, j, 0)[0];
	this->_bufferLeftWall[BUFFER_2D_POS_Y(j)] = flowField.getVelocity().getVector(i+2, j, 0)[1];
}

void VelocityBufferFillStencil::applyRightWall(FlowField & flowField, int i, int j) {
	this->_bufferRightWall[BUFFER_2D_POS_X(j)] = flowField.getVelocity().getVector(i-2, j, 0)[0];
	this->_bufferRightWall[BUFFER_2D_POS_Y(j)] = flowField.getVelocity().getVector(i-1, j, 0)[1];
}

void VelocityBufferFillStencil::applyBottomWall(FlowField & flowField, int i, int j) {
	this->_bufferBottomWall[BUFFER_2D_POS_X(i)] = flowField.getVelocity().getVector(i, j+2, 0)[0];
	this->_bufferBottomWall[BUFFER_2D_POS_Y(i)] = flowField.getVelocity().getVector(i, j+2, 0)[1];
}

void VelocityBufferFillStencil::applyTopWall(FlowField & flowField, int i, int j) {
	this->_bufferTopWall[BUFFER_2D_POS_X(i)] = flowField.getVelocity().getVector(i, j-1, 0)[0];
	this->_bufferTopWall[BUFFER_2D_POS_Y(i)] = flowField.getVelocity().getVector(i, j-2, 0)[1];
}

void VelocityBufferFillStencil::applyLeftWall(FlowField & flowField, int i, int j, int k) {

}

void VelocityBufferFillStencil::applyRightWall(FlowField & flowField, int i, int j, int k) {

}

void VelocityBufferFillStencil::applyBottomWall(FlowField & flowField, int i, int j, int k) {

}

void VelocityBufferFillStencil::applyTopWall(FlowField & flowField, int i, int j, int k) {

}

void VelocityBufferFillStencil::applyFrontWall(FlowField & flowField, int i, int j, int k) {

}

void VelocityBufferFillStencil::applyBackWall(FlowField & flowField, int i, int j, int k) {

}

