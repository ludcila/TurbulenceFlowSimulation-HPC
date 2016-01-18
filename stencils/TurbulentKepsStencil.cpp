#include "StencilFunctions.h"
#include "TurbulentKepsStencil.h"


/* Note: not checking for obstacles! */

TurbulentKepsStencil::TurbulentKepsStencil ( const Parameters & parameters ) : FieldStencil<TurbulentFlowField> ( parameters ) {}

void TurbulentKepsStencil::apply ( TurbulentFlowField & flowField, int i, int j ){

	const int obstacle = flowField.getFlags().getValue(i, j);

	/* fluid cells */
	if ((obstacle & OBSTACLE_SELF) == 0) {

		loadLocalVelocity2D           ( flowField  , _localVelocity          , i, j);
		loadLocalMeshsize2D           ( _parameters, _localMeshsize          , i, j);
		loadLocalTurbulentViscosity2D ( flowField  , _localTurbulentViscosity, i, j);
		loadLocalKineticEnergy2D(flowField,  _localK, i, j);
		loadLocalDissipationRate2D(flowField,  _localeps, i, j);
		loadLocalFmu2D(flowField,  _localfmu, i, j);

		flowField.getKineticEnergyNew().getScalar(i, j) = 
			RHSK2D( _localVelocity, _localMeshsize, _localK,_localeps,  _localTurbulentViscosity, _parameters.timestep.dt, _parameters.turbulence.gamma);
	
		flowField.getDissipationRateNew().getScalar(i, j) =
			RHSeps2D( flowField , _localVelocity, _localMeshsize, _localK, _localeps, _localTurbulentViscosity, _parameters, _localfmu, i,  j, _parameters.timestep.dt, _parameters.turbulence.gamma);

	/* obstacle cells */
	} else if ((obstacle & OBSTACLE_SELF) == 1){
	
		flowField.getKineticEnergyNew().getScalar(i, j) = 0;
	
		/* if top cell is fluid */
		if ((obstacle & OBSTACLE_TOP) == 0){
			flowField.getDissipationRateNew().getScalar(i, j) = flowField.getDissipationRate().getScalar(i, j+1);
		}
		/* if bottom cell is fluid */
		if ((obstacle & OBSTACLE_BOTTOM) == 0){
			flowField.getDissipationRateNew().getScalar(i, j) = flowField.getDissipationRate().getScalar(i, j-1);
		}
		/* if right cell is fluid */
		if ((obstacle & OBSTACLE_RIGHT) == 0){
			flowField.getDissipationRateNew().getScalar(i, j) = flowField.getDissipationRate().getScalar(i+1, j);
		}
		/* if left cell is fluid */
		if ((obstacle & OBSTACLE_LEFT) == 0){
			flowField.getDissipationRateNew().getScalar(i, j) = flowField.getDissipationRate().getScalar(i-1, j);
		}
		/* if left cell is fluid */
		if ((obstacle & OBSTACLE_TOP) == 0 && (obstacle & OBSTACLE_RIGHT) == 0){
			flowField.getDissipationRateNew().getScalar(i, j) = (flowField.getDissipationRate().getScalar(i-1, j) + flowField.getDissipationRate().getScalar(i-1, j)) * 0.5;
		}
		
	}

}



void TurbulentKepsStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ){

}
