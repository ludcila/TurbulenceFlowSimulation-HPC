#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "StencilFunctions.h"
#include "../TurbulentFlowField.h"
#include "TurbulentViscosityKepsStencil.h"
#include "Iterators.h"
#include <iomanip>
#include <mpi.h>



TurbulentViscosityKepsStencil::TurbulentViscosityKepsStencil ( const Parameters & parameters ) :
	FieldStencil<TurbulentFlowField> ( parameters )
{

}

TurbulentViscosityKepsStencil::~TurbulentViscosityKepsStencil () {

}

void TurbulentViscosityKepsStencil::apply ( TurbulentFlowField & flowField, int i, int j ) {

    loadLocalVelocity2D           ( flowField  , _localVelocity          , i, j);
    loadLocalMeshsize2D           ( _parameters, _localMeshsize          , i, j);
    loadLocalTurbulentViscosity2D ( flowField  , _localTurbulentViscosity, i, j);
    loadLocalKineticEnergy2D      ( flowField  ,_localK  , i, j);
    loadLocalDissipationRate2D    ( flowField  ,_localeps, i, j);
    loadLocalFmu2D                ( flowField  ,_localfmu, i, j);

flowField.getTurbulentViscosity().getScalar(i, j) = 1.0/_parameters.flow.Re+	
	  _parameters.turbulence.cmu
  	 *flowField.getKineticEnergy().getScalar(i,j)
	 *flowField.getKineticEnergy().getScalar(i,j)
	 /flowField.getDissipationRate().getScalar(i,j);	

}



void TurbulentViscosityKepsStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ) {


 const int obstacle = flowField.getFlags().getValue(i, j, k);
    

    if ((obstacle & OBSTACLE_SELF) == 0){   // If the cell is fluid

        loadLocalVelocity3D		( flowField , _localVelocity	      , i, j, k);
        loadLocalMeshsize3D		(_parameters, _localMeshsize	      , i, j, k);
        loadLocalTurbulentViscosity3D 	( flowField , _localTurbulentViscosity, i, j, k);
	loadLocalKineticEnergy3D	( flowField , _localK  , i, j, k);
        loadLocalDissipationRate3D	( flowField , _localeps, i, j, k);
    	loadLocalFmu3D                  ( flowField , _localfmu, i, j, k);

        flowField.getTurbulentViscosity().getScalar(i, j, k) =
	 _parameters.turbulence.cmu
	*computefmu3D( flowField , _localK, _localeps, _parameters, i,  j, k)
	*RHSK3D( _localVelocity, _localMeshsize, _localK,_localeps,  _localTurbulentViscosity, _parameters.timestep.dt )
	*RHSK3D( _localVelocity, _localMeshsize, _localK,_localeps,  _localTurbulentViscosity, _parameters.timestep.dt)
	/RHSeps3D( flowField , _localVelocity, _localMeshsize, _localK, _localeps, _localTurbulentViscosity, _parameters, _localfmu , i, j, k,  _parameters.timestep.dt );
    }

}
