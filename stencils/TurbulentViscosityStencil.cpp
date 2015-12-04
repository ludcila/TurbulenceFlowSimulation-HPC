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
FLOAT deltay;
loadLocalVelocity2D(  flowField, _localVelocity, i, j);
loadLocalMeshsize2D(_parameters, _localMeshsize, i, j);
FLOAT re=this->_parameters.flow.Re;



//i set the thickness of the boundary layer for the channel case
if(this->_parameters.simulation.scenario=="channel"){

	if (this->_parameters.turbulence.boundary_layer_equation=="laminar"){
	deltay=4.91*this->_parameters.meshsize->getPosX(i,j)/(pow(re, 1/2.0));
	}

	else 	deltay=0.382*this->_parameters.meshsize->getPosX(i,j)/(pow(re, 1/5.0));

}
 
//if we are in the cavity case, we can apply empirical laws like the Blasius', 'cause they are meant for flate plates, so we just need wall distance
if(this->_parameters.simulation.scenario=="cavity"){
	
	lm=0.41*PetscMin(this->_parameters.meshsize->getPosX(i,j),this->_parameters.meshsize->getPosY(i,j));

	if(i>this->_parameters.geometry.sizeX/2.0 && j<=this->_parameters.geometry.sizeY/2.0)
	lm=0.41*PetscMin(this->_parameters.geometry.lengthX-this->_parameters.meshsize->getPosX(i,j),this->_parameters.meshsize->getPosY(i,j));

	if(i>this->_parameters.geometry.sizeX/2.0 && j>this->_parameters.geometry.sizeY/2.0)
	lm=0.41*PetscMin(this->_parameters.geometry.lengthX-this->_parameters.meshsize->getPosX(i,j),this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j));

	if(i<=this->_parameters.geometry.sizeX/2.0 && j>this->_parameters.geometry.sizeY/2.0)
	lm=0.41*PetscMin(this->_parameters.meshsize->getPosX(i,j),this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j));


}

//now we determine lm for the channel case

if(this->_parameters.simulation.scenario=="channel"){

	if(j>this->_parameters.geometry.sizeY/2.0)
	lm=PetscMin(0.09*deltay,0.41*(this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j)));

	 if(j<=this->_parameters.geometry.sizeY/2.0 && this->_parameters.bfStep.xRatio!=0)
		{
		if(this->_parameters.meshsize->getPosX(i,j)<=this->_parameters.bfStep.xRatio*this->_parameters.geometry.lengthX)
		lm=PetscMin(0.09*deltay,0.41*(-this->_parameters.meshsize->getPosY(i,j)+this->_parameters.bfStep.yRatio*this->_parameters.geometry.lengthY/2.0));

		else lm=0.41*this->_parameters.meshsize->getPosY(i,j);
		}
	
	else lm=PetscMin(0.09*deltay,0.41*this->_parameters.meshsize->getPosY(i,j));
}


flowField.getTurbulentViscosity().getScalar(i,j)=lm*lm*sqrt(2*computeSdotS2D(_localVelocity, _localMeshsize));


}



void TurbulentViscosityStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ) {
FLOAT lm;
FLOAT deltay;
FLOAT deltaz;
loadLocalVelocity3D(  flowField, _localVelocity, i, j , k);
loadLocalMeshsize3D(_parameters, _localMeshsize, i, j , k);
FLOAT re=this->_parameters.flow.Re;



//i set the thickness of the boundary layer for the channel case
if(this->_parameters.simulation.scenario=="channel"){

	if (this->_parameters.turbulence.boundary_layer_equation=="laminar"){
	deltay=4.91*this->_parameters.meshsize->getPosX(i,j,k)/(pow(re, 1/2.0));
	}

	else 	deltay=0.382*this->_parameters.meshsize->getPosX(i,j,k)/(pow(re, 1/5.0));
deltaz=deltay;
}
 

//if we are in the cavity case, we can apply empirical laws like the Blasius', 'cause they are meant for flate plates, so we just need wall distance
if(this->_parameters.simulation.scenario=="cavity"){
	
	if(i<=this->_parameters.geometry.sizeX/2.0 && j<=this->_parameters.geometry.sizeY/2.0 && k<=this->_parameters.geometry.sizeZ/2.0)
	lm=0.41*PetscMin(this->_parameters.meshsize->getPosZ(i,j,k),PetscMin(this->_parameters.meshsize->getPosX(i,j,k),this->_parameters.meshsize->getPosY(i,j,k)));

	if(i>this->_parameters.geometry.sizeX/2.0 && j<=this->_parameters.geometry.sizeY/2.0 && k<=this->_parameters.geometry.sizeZ/2.0)
	lm=0.41*PetscMin(this->_parameters.meshsize->getPosZ(i,j,k),PetscMin(this->_parameters.geometry.lengthX-this->_parameters.meshsize->getPosX(i,j,k),this->_parameters.meshsize->getPosY(i,j,k)));

	if(i>this->_parameters.geometry.sizeX/2.0 && j>this->_parameters.geometry.sizeY/2.0 && k<=this->_parameters.geometry.sizeZ/2.0)
	lm=0.41*PetscMin(this->_parameters.meshsize->getPosZ(i,j,k),PetscMin(this->_parameters.geometry.lengthX-this->_parameters.meshsize->getPosX(i,j,k),this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j,k)));

	if(i<=this->_parameters.geometry.sizeX/2.0 && j>this->_parameters.geometry.sizeY/2.0 && k<=this->_parameters.geometry.sizeZ/2.0)
	lm=0.41*PetscMin(this->_parameters.meshsize->getPosZ(i,j,k),PetscMin(this->_parameters.meshsize->getPosX(i,j,k),this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j,k)));

	if(i<=this->_parameters.geometry.sizeX/2.0 && j<=this->_parameters.geometry.sizeY/2.0 && k>this->_parameters.geometry.sizeZ/2.0)
	lm=0.41*PetscMin(this->_parameters.geometry.lengthZ-this->_parameters.meshsize->getPosZ(i,j,k),PetscMin(this->_parameters.meshsize->getPosX(i,j,k),this->_parameters.meshsize->getPosY(i,j,k)));

	if(i>this->_parameters.geometry.sizeX/2.0 && j<=this->_parameters.geometry.sizeY/2.0 && k>this->_parameters.geometry.sizeZ/2.0)
	lm=0.41*PetscMin(this->_parameters.geometry.lengthZ-this->_parameters.meshsize->getPosZ(i,j,k),PetscMin(this->_parameters.geometry.lengthX-this->_parameters.meshsize->getPosX(i,j,k),this->_parameters.meshsize->getPosY(i,j,k)));

	if(i>this->_parameters.geometry.sizeX/2.0 && j>this->_parameters.geometry.sizeY/2.0 && k>this->_parameters.geometry.sizeZ/2.0)
	lm=0.41*PetscMin(this->_parameters.geometry.lengthZ-this->_parameters.meshsize->getPosZ(i,j,k),PetscMin(this->_parameters.geometry.lengthX-this->_parameters.meshsize->getPosX(i,j,k),this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j,k)));

	if(i<=this->_parameters.geometry.sizeX/2.0 && j>this->_parameters.geometry.sizeY/2.0 && k>this->_parameters.geometry.sizeZ/2.0)
	lm=0.41*PetscMin(this->_parameters.geometry.lengthZ-this->_parameters.meshsize->getPosZ(i,j,k),PetscMin(this->_parameters.meshsize->getPosX(i,j,k),this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j,k)));
}

//now we determine lm for the channel case

if(this->_parameters.simulation.scenario=="channel" && k<=this->_parameters.geometry.sizeZ/2.0){

	if(j>this->_parameters.geometry.sizeY/2.0)
	lm=PetscMin(0.09*deltay,0.41*PetscMin(this->_parameters.meshsize->getPosZ(i,j,k),(this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j,k))));

	 if(j<=this->_parameters.geometry.sizeY/2.0 && this->_parameters.bfStep.xRatio!=0)
		{
		if(this->_parameters.meshsize->getPosX(i,j,k)<=this->_parameters.bfStep.xRatio*this->_parameters.geometry.lengthX)
		lm=PetscMin(0.09*deltay,0.41*PetscMin(this->_parameters.meshsize->getPosZ(i,j,k),(-this->_parameters.meshsize->getPosY(i,j,k)+this->_parameters.bfStep.yRatio*this->_parameters.geometry.lengthY/2.0)));

		else lm=0.41*PetscMin(this->_parameters.meshsize->getPosZ(i,j,k),this->_parameters.meshsize->getPosY(i,j,k));
		}
	
	else lm=PetscMin(0.09*deltay,0.41*PetscMin(this->_parameters.meshsize->getPosZ(i,j,k),this->_parameters.meshsize->getPosY(i,j,k)));
}

if(this->_parameters.simulation.scenario=="channel" && k>this->_parameters.geometry.sizeZ/2.0){

	if(j>this->_parameters.geometry.sizeY/2.0)
	lm=PetscMin(0.09*deltay,0.41*PetscMin(this->_parameters.geometry.lengthZ-this->_parameters.meshsize->getPosZ(i,j,k),(this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j,k))));

	 if(j<=this->_parameters.geometry.sizeY/2.0 && this->_parameters.bfStep.xRatio!=0)
		{
		if(this->_parameters.meshsize->getPosX(i,j,k)<=this->_parameters.bfStep.xRatio*this->_parameters.geometry.lengthX)
		lm=PetscMin(0.09*deltay,0.41*PetscMin(this->_parameters.geometry.lengthZ-this->_parameters.meshsize->getPosZ(i,j,k),(-this->_parameters.meshsize->getPosY(i,j,k)+this->_parameters.bfStep.yRatio*this->_parameters.geometry.lengthY/2.0)));

		else lm=0.41*PetscMin(this->_parameters.geometry.lengthZ-this->_parameters.meshsize->getPosZ(i,j,k),this->_parameters.meshsize->getPosY(i,j,k));
		}
	
	else lm=PetscMin(0.09*deltay,0.41*PetscMin(this->_parameters.geometry.lengthZ-this->_parameters.meshsize->getPosZ(i,j,k),this->_parameters.meshsize->getPosY(i,j,k)));
}


flowField.getTurbulentViscosity().getScalar(i,j,k)=lm*lm*sqrt(2*computeSdotS2D(_localVelocity, _localMeshsize));


}


