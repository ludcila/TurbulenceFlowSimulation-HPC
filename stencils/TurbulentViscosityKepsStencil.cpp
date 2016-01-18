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

	const int obstacle = flowField.getFlags().getValue(i, j);

	if ((obstacle & OBSTACLE_SELF) == 0){   // If the cell is fluid

		flowField.getTurbulentViscosity().getScalar(i, j) = 
			  _parameters.turbulence.cmu
		  	 * flowField.getKineticEnergy().getScalar(i,j)
			 * flowField.getKineticEnergy().getScalar(i,j)
			 / flowField.getDissipationRate().getScalar(i,j);
			 
	 }

}



void TurbulentViscosityKepsStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ) {

	const int obstacle = flowField.getFlags().getValue(i, j, k);

	if ((obstacle & OBSTACLE_SELF) == 0){   // If the cell is fluid
	}

}
