#include "MinNUStencil.h"
#include <algorithm>
#include <math.h>


MinNUStencil::MinNUStencil (const Parameters & parameters) :
    FieldStencil<TurbulentFlowField> (parameters), BoundaryStencil<TurbulentFlowField> (parameters) {
    reset();
}

void MinNUStencil::apply (TurbulentFlowField & flowField, int i, int j){
    cellMinValue(flowField, i, j);
}

void MinNUStencil::apply (TurbulentFlowField & flowField, int i, int j, int k){
    cellMinValue(flowField, i, j, k);
}

void MinNUStencil::applyLeftWall   ( TurbulentFlowField & flowField, int i, int j ){
    cellMinValue(flowField, i, j);
}

void MinNUStencil::applyRightWall  ( TurbulentFlowField & flowField, int i, int j ){
    cellMinValue(flowField, i, j);
}

void MinNUStencil::applyBottomWall ( TurbulentFlowField & flowField, int i, int j ){
    cellMinValue(flowField, i, j);
}

void MinNUStencil::applyTopWall    ( TurbulentFlowField & flowField, int i, int j ){
    cellMinValue(flowField, i, j);
}

void MinNUStencil::applyLeftWall   ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMinValue(flowField, i, j, k);
}

void MinNUStencil::applyRightWall  ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMinValue(flowField, i, j, k);
}

void MinNUStencil::applyBottomWall ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMinValue(flowField, i, j, k);
}

void MinNUStencil::applyTopWall    ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMinValue(flowField, i, j, k);
}

void MinNUStencil::applyFrontWall  ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMinValue(flowField, i, j, k);
}

void MinNUStencil::applyBackWall   ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMinValue(flowField, i, j, k);
}


void MinNUStencil::cellMinValue(TurbulentFlowField & flowField, int i, int j){
    FLOAT * viscosity = flowField.getTurbulentViscosity().getScalar(i, j);
   
    if (viscosity > _minValue){
        _minValue = viscosity;
    }
   
}

void MinNUStencil::cellMinValue(TurbulentFlowField & flowField, int i, int j, int k){
    FLOAT * viscosity= flowField.getTurbulentViscosity().getScalar(i, j);
   
    if (viscosity > _minValue){
        _minValue = viscosity;
    }
}

void MinNUStencil::reset () {
    _minValue = 0;
}

const FLOAT  MinNUStencil::getMinValue() const{
    return _minValue;
}
