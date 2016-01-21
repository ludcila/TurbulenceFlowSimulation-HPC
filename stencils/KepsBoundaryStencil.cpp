#include "KepsBoundaryStencil.h"

KepsBoundaryStencil::KepsBoundaryStencil ( const Parameters & parameters ) :
    BoundaryStencil<TurbulentFlowField> ( parameters ),
    _parameters(parameters) {}

// 2D stencils

void KepsBoundaryStencil::applyLeftWall ( TurbulentFlowField & flowField, int i, int j ){

FLOAT Ret;
FLOAT Redelta;


	const int obstacle = flowField.getFlags().getValue(i, j);
	if((obstacle & OBSTACLE_SELF) == 0) {


	//FLOAT L=0.07*_parameters.geometry.lengthY;
	//FLOAT I= 0.16 *pow(_parameters.flow.Re,-1/8);
	//FLOAT I=sqrt(2/3*flowField.getKineticEnergy().getScalar(i+1, j));
	//FLOAT beta=flowField.getTurbulentViscosity().getScalar(i+1, j)*_parameters.flow.Re;
	
	//flowField.getKineticEnergy().getScalar(i, j)=1.5*pow(flowField.getVelocity().getVector(i,j)[0]*I,2);
	//flowField.getDissipationRate().getScalar(i, j) =0.09*pow(flowField.getKineticEnergy().getScalar(i, j),2)/(1/_parameters.flow.Re*beta);
    //flowField.getKineticEnergy().getScalar(i, j) = 0.003 * _parameters.walls.vectorLeft[0];
     //flowField.getKineticEnergy().getScalar(i, j) = 0.003*pow(flowField.getVelocity().getVector(i,j)[0],2);
    //flowField.getDissipationRate().getScalar(i, j) = _parameters.turbulence.cmu*pow( flowField.getKineticEnergy().getScalar(i, j), 1.5 )/(0.03*_parameters.geometry.lengthY);

	flowField.getKineticEnergy().getScalar(i, j)=1.5* pow(_parameters.walls.vectorLeft[0]*0.05,2);
	flowField.getDissipationRate().getScalar(i, j) =0.09/(0.038*_parameters.geometry.lengthY*(1-_parameters.bfStep.yRatio))*pow(flowField.getKineticEnergy().getScalar(i, j),1.5);


 Ret= flowField.getKineticEnergy().getScalar(i, j)*flowField.getKineticEnergy().getScalar(i, j)*_parameters.flow.Re/flowField.getDissipationRate().getScalar(i, j);

 Redelta=sqrt(flowField.getKineticEnergy().getScalar(i, j))*flowField.getWallDistance(i, j)*_parameters.flow.Re;
   	flowField.getTurbulentViscosity().getScalar(i, j) = 
    		_parameters.turbulence.cmu 
		*(1-exp(-0.0165*Redelta))*(1-exp(-0.0165*Redelta))*(1+20.5/Ret) 
    		* flowField.getKineticEnergy().getScalar(i, j)
    		* flowField.getKineticEnergy().getScalar(i, j)
    		/ flowField.getDissipationRate().getScalar(i, j);
	}
}


void KepsBoundaryStencil::applyRightWall ( TurbulentFlowField & flowField, int i, int j ){
	flowField.getKineticEnergy().getScalar(i, j) =  flowField.getKineticEnergy().getScalar(i-1, j);
	flowField.getDissipationRate().getScalar(i, j) =  flowField.getDissipationRate().getScalar(i-1, j);
	flowField.getTurbulentViscosity().getScalar(i, j) = flowField.getTurbulentViscosity().getScalar(i-1, j);
}


void KepsBoundaryStencil::applyBottomWall ( TurbulentFlowField & flowField, int i, int j ){

	const int obstacle = flowField.getFlags().getValue(i, j);
	if((obstacle & OBSTACLE_SELF) == 0) {
		flowField.getKineticEnergy().getScalar(i, j) = -flowField.getKineticEnergy().getScalar(i, j+1) * _parameters.meshsize->getDy(i, j)/_parameters.meshsize->getDy(i, j+1) ;
    	flowField.getDissipationRate().getScalar(i, j) =  flowField.getDissipationRate().getScalar(i, j+1);
    	flowField.getTurbulentViscosity().getScalar(i, j) = -flowField.getTurbulentViscosity().getScalar(i, j+1) * _parameters.meshsize->getDy(i, j)/_parameters.meshsize->getDy(i, j+1) ;
	}
	
}

void KepsBoundaryStencil::applyTopWall ( TurbulentFlowField & flowField, int i, int j ){
	flowField.getKineticEnergy().getScalar(i, j) = -flowField.getKineticEnergy().getScalar(i, j-1) * _parameters.meshsize->getDy(i, j)/_parameters.meshsize->getDy(i, j-1);
	flowField.getDissipationRate().getScalar(i, j) =  flowField.getDissipationRate().getScalar(i, j-1);
	flowField.getTurbulentViscosity().getScalar(i, j) = -flowField.getTurbulentViscosity().getScalar(i, j-1) * _parameters.meshsize->getDy(i, j)/_parameters.meshsize->getDy(i, j-1);
}


// 3D stencils

void KepsBoundaryStencil::applyLeftWall ( TurbulentFlowField & flowField, int i, int j, int k ){


	flowField.getKineticEnergy().getScalar(i, j, k)=1.5* pow(_parameters.walls.vectorLeft[0]*0.05,2);
	flowField.getDissipationRate().getScalar(i, j, k) =0.09/(0.038*_parameters.geometry.lengthY*(1-_parameters.bfStep.yRatio)*_parameters.geometry.lengthZ*2/(_parameters.geometry.lengthY*(1-_parameters.bfStep.yRatio)+_parameters.geometry.lengthZ))*pow(flowField.getKineticEnergy().getScalar(i, j, k),1.5);
    	flowField.getTurbulentViscosity().getScalar(i, j, k) = 
    		_parameters.turbulence.cmu 
    		* flowField.getKineticEnergy().getScalar(i, j, k)
    		* flowField.getKineticEnergy().getScalar(i, j, k)
    		/ flowField.getDissipationRate().getScalar(i, j, k);
}


void KepsBoundaryStencil::applyRightWall ( TurbulentFlowField & flowField, int i, int j , int k ){
     	flowField.getKineticEnergy().getScalar(i, j,k) =  flowField.getKineticEnergy().getScalar(i-1, j,k);
     	flowField.getDissipationRate().getScalar(i, j,k) =  flowField.getDissipationRate().getScalar(i-1, j,k);
     	flowField.getTurbulentViscosity().getScalar(i, j,k) = flowField.getTurbulentViscosity().getScalar(i-1, j,k);
}


void KepsBoundaryStencil::applyBottomWall ( TurbulentFlowField & flowField, int i, int j, int k ){
	flowField.getKineticEnergy().getScalar(i, j, k) =-flowField.getKineticEnergy().getScalar(i, j+1, k) * _parameters.meshsize->getDy(i, j, k)/_parameters.meshsize->getDy(i, j+1, k) ;
    	flowField.getDissipationRate().getScalar(i, j, k) =  flowField.getDissipationRate().getScalar(i, j+1, k);
  	flowField.getTurbulentViscosity().getScalar(i, j,k ) = -flowField.getTurbulentViscosity().getScalar(i, j+1, k) * _parameters.meshsize->getDy(i, j, k)/_parameters.meshsize->getDy(i, j+1, k) ;
    		
}


void KepsBoundaryStencil::applyTopWall ( TurbulentFlowField & flowField, int i, int j, int k ){

	flowField.getKineticEnergy().getScalar(i, j, k) =-flowField.getKineticEnergy().getScalar(i, j-1, k) * _parameters.meshsize->getDy(i, j, k)/_parameters.meshsize->getDy(i, j-1, k) ;
    	flowField.getDissipationRate().getScalar(i, j, k) =  flowField.getDissipationRate().getScalar(i, j-1, k);
    	flowField.getTurbulentViscosity().getScalar(i, j, k) =  -flowField.getTurbulentViscosity().getScalar(i, j-1, k) * _parameters.meshsize->getDy(i, j,k)/_parameters.meshsize->getDy(i, j-1,k );
}


void KepsBoundaryStencil::applyFrontWall ( TurbulentFlowField & flowField, int i, int j, int k ){

	flowField.getKineticEnergy().getScalar(i, j, k) =-flowField.getKineticEnergy().getScalar(i, j, k+1) * _parameters.meshsize->getDz(i, j, k)/_parameters.meshsize->getDz(i, j, k+1); 
    	flowField.getDissipationRate().getScalar(i, j, k) =  flowField.getDissipationRate().getScalar(i, j, k+1);
  	flowField.getTurbulentViscosity().getScalar(i, j, k) =  -flowField.getTurbulentViscosity().getScalar(i, j,k+1) * _parameters.meshsize->getDz(i, j,k)/_parameters.meshsize->getDz(i, j,k+1);
}


void KepsBoundaryStencil::applyBackWall ( TurbulentFlowField & flowField, int i, int j, int k ){

	flowField.getKineticEnergy().getScalar(i, j, k) =-flowField.getKineticEnergy().getScalar(i, j, k-1) * _parameters.meshsize->getDz(i, j, k)/_parameters.meshsize->getDz(i, j, k-1); 
    	flowField.getDissipationRate().getScalar(i, j, k) =  flowField.getDissipationRate().getScalar(i, j, k-1);
 	flowField.getTurbulentViscosity().getScalar(i, j, k) =  -flowField.getTurbulentViscosity().getScalar(i, j,k-1) * _parameters.meshsize->getDz(i, j,k)/_parameters.meshsize->getDz(i, j,k-1);

}
