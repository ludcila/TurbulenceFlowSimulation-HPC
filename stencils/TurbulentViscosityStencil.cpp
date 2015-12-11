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
FLOAT xpos = this->_parameters.meshsize->getPosX(i,j) + this->_parameters.meshsize->getDx(i,j)/2;
int sizeY = this->_parameters.geometry.sizeY;
//FLOAT re=flowField.getVelocity().getVector(i,sizeY/2)[0]*xpos*this->_parameters.flow.Re;
FLOAT re=flowField.getCenterLineVelocity()[i]*xpos*this->_parameters.flow.Re;
//FLOAT re=flowField.getCenterLineVelocity()[0]*xpos*this->_parameters.flow.Re;
//std::cout << std::endl <<_parameters.parallel.indices[1] << std::endl;



//i set the thickness of the boundary layer for the channel case
if(this->_parameters.simulation.scenario=="channel"){

	if (this->_parameters.turbulence.boundary_layer_equation=="laminar"){
	deltay=4.91*xpos/(pow(re, 1/2.0));
	}

	else 	deltay=0.382*xpos/(pow(re, 1/5.0));

}
 
//if we are in the cavity case, we can apply empirical laws like the Blasius', 'cause they are meant for flate plates, so we just need wall distance
if(this->_parameters.simulation.scenario=="cavity"){
	lm=0.41*flowField.getWallDistance().getScalar(i, j);


}

//now we determine lm for the channel case

if(this->_parameters.simulation.scenario=="channel"){

		FLOAT a = 0.09*deltay;
		FLOAT b = 0.41*flowField.getWallDistance().getScalar(i, j);

		/*
		if(this->_parameters.bfStep.xRatio!=0 && this->_parameters.bfStep.yRatio && xpos>this->_parameters.bfStep.xRatio*this->_parameters.geometry.lengthX && this->_parameters.meshsize->getPosY(i,j)<=this->_parameters.geometry.lengthY/2)
		lm=b;
		
	
		else
		*/
		lm = PetscMin(a, b);
	
}


flowField.getTurbulentViscosity().getScalar(i,j)=lm*lm*sqrt(2*computeSdotS2D(_localVelocity, _localMeshsize));


}



void TurbulentViscosityStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ) {
FLOAT lm;
FLOAT deltay;
FLOAT deltaz;
loadLocalVelocity3D(  flowField, _localVelocity, i, j , k);
loadLocalMeshsize3D(_parameters, _localMeshsize, i, j , k);
FLOAT xpos = this->_parameters.meshsize->getPosX(i,j,k) + this->_parameters.meshsize->getDx(i,j,k);
int sizeY = this->_parameters.geometry.sizeY;
int sizeZ = this->_parameters.geometry.sizeZ;
FLOAT re=flowField.getVelocity().getVector(i,sizeY/2,sizeZ/2)[0]*xpos*this->_parameters.flow.Re;
//FLOAT re=this->_parameters.flow.Re;



//i set the thickness of the boundary layer for the channel case
if(this->_parameters.simulation.scenario=="channel"){

	if (this->_parameters.turbulence.boundary_layer_equation=="laminar"){
	deltay=4.91*xpos/(pow(re, 1/2.0));
	}

	else 	deltay=0.382*xpos/(pow(re, 1/5.0));
deltaz=deltay;
}
 

//if we are in the cavity case, we can apply empirical laws like the Blasius', 'cause they are meant for flate plates, so we just need wall distance
if(this->_parameters.simulation.scenario=="cavity"){
	
	lm=0.41*flowField.getWallDistance().getScalar(i, j,k);
}

//now we determine lm for the channel case

if(this->_parameters.simulation.scenario=="channel"){

    		FLOAT a = 0.09*deltay;
		FLOAT b = 0.41*flowField.getWallDistance().getScalar(i, j, k);

		if(this->_parameters.bfStep.xRatio!=0 && this->_parameters.bfStep.yRatio && xpos>this->_parameters.bfStep.xRatio*this->_parameters.geometry.lengthX && this->_parameters.meshsize->getPosY(i,j,k)<=this->_parameters.geometry.lengthY/2)
		lm=b;
		
		else
		lm = PetscMin(a, b);
	


}


//std::cout << i << " " << j << " " << k << " " << ss << std::endl;
flowField.getTurbulentViscosity().getScalar(i,j,k)=lm*lm*sqrt(2*computeSdotS3D(_localVelocity, _localMeshsize));


}


