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
const int obstacleleft_down = flowField.getFlags().getValue(this->_parameters.parallel.firstCorner[0],this->_parameters.parallel.firstCorner[1]);
/*const int obstacleright_up = flowField.getFlags().getValue(this->_parameters.parallel.firstcorner(0)+this->_parameters.geometry.sizeX-1,this->_parameters.parallel.firstcorner(1)+this->_parameters.geometry.sizeY-1);*/
FLOAT mindistance;
FLOAT delta;
FLOAT re=this->_parameters.flow.Re;

 
if(this->_parameters.turbulence.boundary_layer_equation=="laminar")
{
	if((obstacleleft_down & OBSTACLE_LEFT)==1)  delta=1000;

	else delta=4.91*this->_parameters.meshsize->getPosX(i,j)/(pow(re, 1/2));
}

else {
	if((obstacleleft_down & OBSTACLE_LEFT)==1)  delta=1000;

	else delta=0.382*this->_parameters.meshsize->getPosX(i,j)/(pow(re, 1/5));
}
	

if((obstacleleft_down & OBSTACLE_LEFT)==1)
	{
	if(j<=i || j<this->_parameters.geometry.sizeX-i)
		if( j<=this->_parameters.geometry.sizeY/2) mindistance=this->_parameters.meshsize->getPosY(i,j);
	
	else if(this->_parameters.geometry.sizeY-j<=i || this->_parameters.geometry.sizeY-j<=this->_parameters.geometry.sizeX-i)
		if( j>this->_parameters.geometry.sizeY/2) mindistance=this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j);

	else{
		if(i<=this->_parameters.geometry.sizeX/2)  mindistance=this->_parameters.meshsize->getPosX(i,j);

		else mindistance=this->_parameters.geometry.lengthX-this->_parameters.meshsize->getPosX(i,j);

}

else {

	if(this->_parameters.meshsize->getPosX(i,j)<this->_parameters.bfStep.xRatio){
		if(this->_parameters.meshsize->getPosY(i,j)-this->_parameters.bfStep.yRatio<=this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j))   mindistance=this->_parameters.meshsize->getPosY(i,j)-this->_parameters.bfStep.yRatio;
		else mindistance=this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j);
		
	}

	else{
		if(j<=this->_parameters.geometry.sizeY/2) mindistance=this->_parameters.meshsize->getPosY(i,j);
		
		else mindistance=this->_parameters.geometry.lengthY-this->_parameters.meshsize->getPosY(i,j);
	}

}



if(0.41*this->_parameters.meshsize->getPosY(i,j) < 0.09*4.91* this->_parameters.meshsize->getPosX(i,j)/(sqrt(this->_parameters.flow.Re)))
lm=0.41*this->_parameters.meshsize->getPosY(i,j);

else 
lm= 0.09*4.91* this->_parameters.meshsize->getPosX(i,j)/(sqrt(this->_parameters.flow.Re));

flowField.getTurbulentViscosity().getScalar(i,j)=lm*lm*sqrt(2*computeSdotS2D(_localVelocity, _localMeshsize));

}

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


