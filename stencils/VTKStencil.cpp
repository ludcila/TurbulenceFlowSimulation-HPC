#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include "VTKStencil.h"
#include <string>
#include <fstream>
#include <sstream>

VTKStencil::VTKStencil ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) {

}

void VTKStencil::apply ( FlowField & flowField, int i, int j ) {

}

void VTKStencil::apply ( FlowField & flowField, int i, int j, int k ) {

}

void VTKStencil::write ( FlowField & flowField, int timeStep ) {

}
