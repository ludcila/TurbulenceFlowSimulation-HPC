#include "MaxNuStencil.h"
#include <algorithm>
#include <math.h>


MaxNuStencil::MaxNuStencil (const Parameters & parameters) :
    FieldStencil<TurbulentFlowField> (parameters), BoundaryStencil<TurbulentFlowField> (parameters),_parameters(parameters) {
    reset();
}

void MaxNuStencil::apply (TurbulentFlowField & flowField, int i, int j){
    cellMaxValue(flowField, i, j);
}

void MaxNuStencil::apply (TurbulentFlowField & flowField, int i, int j, int k){
    cellMaxValue(flowField, i, j, k);
}

void MaxNuStencil::applyLeftWall   ( TurbulentFlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxNuStencil::applyRightWall  ( TurbulentFlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxNuStencil::applyBottomWall ( TurbulentFlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxNuStencil::applyTopWall    ( TurbulentFlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxNuStencil::applyLeftWall   ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxNuStencil::applyRightWall  ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxNuStencil::applyBottomWall ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxNuStencil::applyTopWall    ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxNuStencil::applyFrontWall  ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxNuStencil::applyBackWall   ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}


void MaxNuStencil::cellMaxValue(TurbulentFlowField & flowField, int i, int j){
    FLOAT  viscosity = flowField.getTurbulentViscosity().getScalar(i, j) + 1.0 / _parameters.flow.Re;
    FLOAT dx2 = _parameters.meshsize->getDx(i, j) * _parameters.meshsize->getDx(i, j);
    FLOAT dy2 = _parameters.meshsize->getDy(i, j) * _parameters.meshsize->getDy(i, j);
    FLOAT newVal = (2 * viscosity * (1.0 / dx2 + 1.0 / dy2));
    if (newVal > _maxValue){
        _maxValue = newVal;
    }
   
}

void MaxNuStencil::cellMaxValue(TurbulentFlowField & flowField, int i, int j, int k){
    FLOAT  viscosity = flowField.getTurbulentViscosity().getScalar(i, j, k) + 1.0 / _parameters.flow.Re;
    FLOAT dx2 = _parameters.meshsize->getDx(i, j, k) * _parameters.meshsize->getDx(i, j, k);
    FLOAT dy2 = _parameters.meshsize->getDy(i, j, k) * _parameters.meshsize->getDy(i, j, k);
    FLOAT dz2 = _parameters.meshsize->getDz(i, j, k) * _parameters.meshsize->getDz(i, j, k);
    FLOAT newVal = (2 * viscosity * (1.0 / dx2 + 1.0 / dy2 + 1.0 / dz2));
    if (newVal > _maxValue){
        _maxValue = newVal;
    }
}

void MaxNuStencil::reset () {
    _maxValue = 0;
}

FLOAT  MaxNuStencil::getMaxValue() const{
    return _maxValue;
}
