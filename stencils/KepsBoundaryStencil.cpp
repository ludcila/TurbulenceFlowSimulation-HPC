#include "KepsBoundaryStencil.h"

KepsBoundaryStencil::KepsBoundaryStencil ( const Parameters & parameters ) :
    BoundaryStencil<TurbulentFlowField> ( parameters ),
    _parameters(parameters) {}

// 2D stencils

void KepsBoundaryStencil::applyLeftWall ( TurbulentFlowField & flowField, int i, int j ){
    flowField.getKineticEnergy().getScalar(i, j) = 0.003*pow(flowField.getVelocity().getVector(i,j)[0],2);
    flowField.getDissipationRate().getScalar(i, j) = _parameters.turbulence.cmu*pow( flowField.getKineticEnergy().getScalar(i, j),3/2 )/(0.03*_parameters.geometry.lengthY);
}


void KepsBoundaryStencil::applyRightWall ( TurbulentFlowField & flowField, int i, int j ){
     flowField.getKineticEnergy().getScalar(i, j) =  flowField.getKineticEnergy().getScalar(i-1, j);
    flowField.getDissipationRate().getScalar(i, j) =  flowField.getDissipationRate().getScalar(i-1, j);
}


void KepsBoundaryStencil::applyBottomWall ( TurbulentFlowField & flowField, int i, int j ){

 flowField.getKineticEnergy().getScalar(i, j) =-flowField.getKineticEnergy().getScalar(i, j+1) * _parameters.meshsize->getDy(i, j)/_parameters.meshsize->getDy(i, j+1); ;
    flowField.getDissipationRate().getScalar(i, j) =  flowField.getDissipationRate().getScalar(i, j+1);
}

void KepsBoundaryStencil::applyTopWall ( TurbulentFlowField & flowField, int i, int j ){

flowField.getKineticEnergy().getScalar(i, j) =-flowField.getKineticEnergy().getScalar(i, j-1) * _parameters.meshsize->getDy(i, j)/_parameters.meshsize->getDy(i, j-1); ;
    flowField.getDissipationRate().getScalar(i, j) =  flowField.getDissipationRate().getScalar(i, j-1);
}


// 3D stencils

void KepsBoundaryStencil::applyLeftWall ( TurbulentFlowField & flowField, int i, int j, int k ){
 flowField.getKineticEnergy().getScalar(i, j,k) = 0.003*pow(flowField.getVelocity().getVector(i,j,k)[0],2);
    flowField.getDissipationRate().getScalar(i, j,k) = _parameters.turbulence.cmu*pow( flowField.getKineticEnergy().getScalar(i, j,k),3/2 )/(0.03*_parameters.geometry.lengthY*_parameters.geometry.lengthZ*2/(_parameters.geometry.lengthY+_parameters.geometry.lengthZ));
}


void KepsBoundaryStencil::applyRightWall ( TurbulentFlowField & flowField, int i, int j , int k ){
     flowField.getKineticEnergy().getScalar(i, j,k) =  flowField.getKineticEnergy().getScalar(i-1, j,k);
    flowField.getDissipationRate().getScalar(i, j,k) =  flowField.getDissipationRate().getScalar(i-1, j,k);

}


void KepsBoundaryStencil::applyBottomWall ( TurbulentFlowField & flowField, int i, int j, int k ){
flowField.getKineticEnergy().getScalar(i, j, k) =-flowField.getKineticEnergy().getScalar(i, j+1, k) * _parameters.meshsize->getDy(i, j, k)/_parameters.meshsize->getDy(i, j+1, k) ;
    flowField.getDissipationRate().getScalar(i, j, k) =  flowField.getDissipationRate().getScalar(i, j+1, k);
}


void KepsBoundaryStencil::applyTopWall ( TurbulentFlowField & flowField, int i, int j, int k ){

flowField.getKineticEnergy().getScalar(i, j, k) =-flowField.getKineticEnergy().getScalar(i, j-1, k) * _parameters.meshsize->getDy(i, j, k)/_parameters.meshsize->getDy(i, j-1, k); ;
    flowField.getDissipationRate().getScalar(i, j, k) =  flowField.getDissipationRate().getScalar(i, j-1, k);
}


void KepsBoundaryStencil::applyFrontWall ( TurbulentFlowField & flowField, int i, int j, int k ){

flowField.getKineticEnergy().getScalar(i, j, k) =-flowField.getKineticEnergy().getScalar(i, j, k-1) * _parameters.meshsize->getDz(i, j, k)/_parameters.meshsize->getDz(i, j, k-1); ;
    flowField.getDissipationRate().getScalar(i, j, k) =  flowField.getDissipationRate().getScalar(i, j, k-1);

}


void KepsBoundaryStencil::applyBackWall ( TurbulentFlowField & flowField, int i, int j, int k ){

flowField.getKineticEnergy().getScalar(i, j, k) =-flowField.getKineticEnergy().getScalar(i, j, k+1) * _parameters.meshsize->getDz(i, j, k)/_parameters.meshsize->getDz(i, j, k+1); ;
    flowField.getDissipationRate().getScalar(i, j, k) =  flowField.getDissipationRate().getScalar(i, j, k+1);
}
