#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
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


if( this->_parameters.meshsize->getPosY(i,j) < 4.91* this->_parameters.meshsize->getPosX(i,j)/(sqrt(this->_parameters.flow.Re)))
flowField.getTurbulentViscosity().getScalar(i,j)=this->_parameters.meshsize->getPosY(i,j);

else 


flowField.getTurbulentViscosity().getScalar(i,j)=4.91* this->_parameters.meshsize->getPosX(i,j)/(sqrt(this->_parameters.flow.Re));
}





void TurbulentViscosityStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ) {
FLOAT closestwall;

if (this->_parameters.meshsize->getPosY(i,j,k)<this->_parameters.meshsize->getPosZ(i,j,k))
closestwall=this->_parameters.meshsize->getPosY(i,j,k);

else 
closestwall=this->_parameters.meshsize->getPosZ(i,j,k);


if( closestwall < 4.91*this->_parameters.meshsize->getPosX(i,j,k)/(sqrt(this->_parameters.flow.Re)))
flowField.getTurbulentViscosity().getScalar(i,j,k)=closestwall;

else 
flowField.getTurbulentViscosity().getScalar(i,j,k)=4.91* this->_parameters.meshsize->getPosX(i,j,k)/(sqrt(this->_parameters.flow.Re));

}


