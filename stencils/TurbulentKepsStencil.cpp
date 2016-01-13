#include "StencilFunctions.h"
#include "TurbulentKepsStencil.h"


/* Note: not checking for obstacles! */

TurbulentKepsStencil::TurbulentKepsStencil ( const Parameters & parameters ) : FieldStencil<TurbulentFlowField> ( parameters ) {}

void TurbulentKepsStencil::apply ( TurbulentFlowField & flowField, int i, int j ){

    loadLocalVelocity2D           ( flowField  , _localVelocity          , i, j);
    loadLocalMeshsize2D           ( _parameters, _localMeshsize          , i, j);
    loadLocalTurbulentViscosity2D ( flowField  , _localTurbulentViscosity, i, j);
    loadLocalKineticEnergy2D(flowField,  _localK, i, j);
    loadLocalDissipationRate2D(flowField,  _localeps, i, j);

	flowField.getKineticEnergy().getScalar(i, j) = RHSK2D( _localVelocity, _localMeshsize, _localK,_localeps,  _localTurbulentViscosity, _parameters.timestep.dt);
	flowField.getDissipationRate().getScalar(i, j) = RHSeps2D( flowField , _localVelocity, _localMeshsize, _localK, _localeps, _localTurbulentViscosity, _parameters,  i,  j, _parameters.timestep.dt);

	if(i==2 && j==2) {
	    std::cout << flowField.getKineticEnergy().getScalar(i, j)  << "\t" << flowField.getDissipationRate().getScalar(i, j) << std::endl;
	}

}



void TurbulentKepsStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ){

}
