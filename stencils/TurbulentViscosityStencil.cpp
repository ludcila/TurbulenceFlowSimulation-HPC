#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "StencilFunctions.h"
#include "../TurbulentFlowField.h"
#include "TurbulentViscosityStencil.h"
#include "Iterators.h"
#include <iomanip>
#include <mpi.h>

TurbulentViscosityStencil::TurbulentViscosityStencil ( const Parameters & parameters ) : FieldStencil<TurbulentFlowField> ( parameters ) {

}

TurbulentViscosityStencil::~TurbulentViscosityStencil () {
	
}

void TurbulentViscosityStencil::apply ( TurbulentFlowField & flowField, int i, int j ) {
FLOAT lm;
loadLocalVelocity2D(  flowField, _localVelocity, i, j);
loadLocalMeshsize2D(_parameters, _localMeshsize, i, j);

if(0.41*this->_parameters.meshsize->getPosY(i,j) < 0.09*4.91* this->_parameters.meshsize->getPosX(i,j)/(sqrt(this->_parameters.flow.Re)))
lm=0.41*this->_parameters.meshsize->getPosY(i,j);

else 
lm= 0.09*4.91* this->_parameters.meshsize->getPosX(i,j)/(sqrt(this->_parameters.flow.Re));

flowField.getTurbulentViscosity().getScalar(i,j)=lm*lm*sqrt(2*computeSdotS2D(_localVelocity, _localMeshsize));

}





void TurbulentViscosityStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ) {
FLOAT closestwall;
FLOAT lm;
loadLocalVelocity3D(  flowField, _localVelocity, i, j, k);
loadLocalMeshsize3D(_parameters, _localMeshsize, i, j, k);

if (this->_parameters.meshsize->getPosY(i,j,k)<this->_parameters.meshsize->getPosZ(i,j,k))
closestwall=this->_parameters.meshsize->getPosY(i,j,k);

else 
closestwall=this->_parameters.meshsize->getPosZ(i,j,k);


if( 0.41*closestwall < 0.09*4.91*this->_parameters.meshsize->getPosX(i,j,k)/(sqrt(this->_parameters.flow.Re)))
lm=0.41*closestwall;

else 
lm=0.09*4.91* this->_parameters.meshsize->getPosX(i,j,k)/(sqrt(this->_parameters.flow.Re));

flowField.getTurbulentViscosity().getScalar(i,j,k)=lm*lm*sqrt(2*computeSdotS3D(_localVelocity, _localMeshsize));
}


