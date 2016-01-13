#include "FmuStencil.h"

FmuStencil::FmuStencil ( const Parameters & parameters ) :
FieldStencil<TurbulentFlowField> ( parameters ){}


void FmuStencil::apply ( TurbulentFlowField & flowField,  int i, int j ){

    loadLocalMeshsize2D       ( _parameters, _localMeshsize   , i, j);
    loadLocalKineticEnergy2D  ( flowField  ,  _localK         , i, j);
    loadLocalDissipationRate2D( flowField  ,  _localeps       , i, j);


    flowField.getFmu().getScalar(i,j) = computefmu2D(flowField, _localK,_localeps,_parameters,i,j);


}


void FmuStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ){

    loadLocalMeshsize3D       ( _parameters, _localMeshsize  , i, j, k);
    loadLocalKineticEnergy3D  ( flowField  , _localK         , i, j, k);
    loadLocalDissipationRate3D( flowField  , _localeps       , i, j, k);


    flowField.getFmu().getScalar(i,j,k) = computefmu3D(flowField, _localK,_localeps,_parameters,i,j,k);
}
