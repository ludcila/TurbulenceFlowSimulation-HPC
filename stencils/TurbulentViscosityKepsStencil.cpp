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


}

void TurbulentViscosityKepsStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ) {

}
