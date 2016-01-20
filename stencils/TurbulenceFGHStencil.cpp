#include "TurbulenceFGHStencil.h"

TurbulenceFGHStencil::TurbulenceFGHStencil ( const Parameters & parameters ) :
FieldStencil<TurbulentFlowField> ( parameters ){}


void TurbulenceFGHStencil::apply ( TurbulentFlowField & flowField,  int i, int j ){

    const int obstacle = flowField.getFlags().getValue(i, j);

    if ((obstacle & OBSTACLE_SELF) == 0){   // If the cell is fluid
    
		// Load local variables into the local array
		loadLocalVelocity2D           ( flowField  , _localVelocity          , i, j);
		loadLocalMeshsize2D           ( _parameters, _localMeshsize          , i, j);
		loadLocalTurbulentViscosity2D ( flowField  , _localTurbulentViscosity, i, j);
		loadLocalKineticEnergy2D      ( flowField  , _localK 		     ,  i, j);

		FLOAT* const values = flowField.getFGH().getVector(i,j);
		                                  
        if ((obstacle & OBSTACLE_RIGHT) == 0) { // If the right cell is fluid
			//values [0] = computeF2DTurbulence(_localVelocity, _localMeshsize, _localTurbulentViscosity, _parameters, _parameters.timestep.dt);
			values [0] = computeF2DTurbulenceKeps(_localVelocity, _localK, _localMeshsize, _localTurbulentViscosity, _parameters, _parameters.timestep.dt);
        }
        if ((obstacle & OBSTACLE_TOP) == 0) {
			//values [1] = computeG2DTurbulence(_localVelocity, _localMeshsize, _localTurbulentViscosity, _parameters, _parameters.timestep.dt);
			values [1] = computeG2DTurbulenceKeps(_localVelocity, _localK, _localMeshsize, _localTurbulentViscosity, _parameters, _parameters.timestep.dt);
        }
		                                  
	}
}


void TurbulenceFGHStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ){

    const int obstacle = flowField.getFlags().getValue(i, j, k);

    if ((obstacle & OBSTACLE_SELF) == 0){   // If the cell is fluid
	FLOAT * const values = flowField.getFGH().getVector(i,j,k);

        loadLocalVelocity3D		(  flowField  , _localVelocity,           i, j, k);
        loadLocalMeshsize3D		(  _parameters, _localMeshsize,           i, j, k);
        loadLocalTurbulentViscosity3D   (  flowField  , _localTurbulentViscosity, i, j, k);
	loadLocalKineticEnergy3D	(  flowField  , _localK, 		  i, j, k);

        if ((obstacle & OBSTACLE_RIGHT) == 0) { // If the right cell is fluid
            values [0] = computeF3DTurbulenceKeps(_localVelocity, _localK, _localMeshsize, _localTurbulentViscosity, _parameters, _parameters.timestep.dt);
        }
        if ((obstacle & OBSTACLE_TOP) == 0) {
            values [1] = computeG3DTurbulenceKeps(_localVelocity, _localK, _localMeshsize, _localTurbulentViscosity, _parameters, _parameters.timestep.dt);
        }
        if ((obstacle & OBSTACLE_BACK) == 0) {
            values [2] = computeH3DTurbulenceKeps(_localVelocity, _localK, _localMeshsize, _localTurbulentViscosity, _parameters, _parameters.timestep.dt);
        }
    }
}
