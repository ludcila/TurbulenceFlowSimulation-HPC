#include "StencilFunctions.h"
#include "TurbulentKepsStencil.h"


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

	}

}



void TurbulentKepsStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ){

	const int obstacle = flowField.getFlags().getValue(i, j, k);

	/* fluid cells */
	if ((obstacle & OBSTACLE_SELF) == 0) {

		loadLocalVelocity3D           ( flowField  , _localVelocity          , i, j, k);
		loadLocalMeshsize3D           ( _parameters, _localMeshsize          , i, j,k );
		loadLocalTurbulentViscosity3D ( flowField  , _localTurbulentViscosity, i, j, k);
		loadLocalKineticEnergy3D(flowField,  _localK, i, j ,k);
		loadLocalDissipationRate3D(flowField,  _localeps, i, j, k);
		loadLocalFmu3D(flowField,  _localfmu, i, j, k);

		flowField.getKineticEnergyNew().getScalar(i, j, k) = 
			RHSK3D( _localVelocity, _localMeshsize, _localK,_localeps,  _localTurbulentViscosity, _parameters.timestep.dt, _parameters.turbulence.gamma);
	
		flowField.getDissipationRateNew().getScalar(i, j, k) =
			RHSeps3D( flowField , _localVelocity, _localMeshsize, _localK, _localeps, _localTurbulentViscosity, _parameters, _localfmu, i,  j, k, _parameters.timestep.dt, _parameters.turbulence.gamma);

	}
}

