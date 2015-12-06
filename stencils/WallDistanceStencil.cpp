#include <algorithm>
#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../TurbulentFlowField.h"
#include "VTKStencil.h"
#include "Iterators.h"
#include "WallDistanceStencil.h"

WallDistanceStencil::WallDistanceStencil ( const Parameters & parameters ) : 
	FieldStencil<TurbulentFlowField> ( parameters ),
	_lengthX(parameters.geometry.lengthX),
	_lengthY(parameters.geometry.lengthY),
	_lengthZ(parameters.geometry.lengthZ)
{
}

WallDistanceStencil::~WallDistanceStencil () {
}

void WallDistanceStencil::apply ( TurbulentFlowField & flowField, int i, int j ) {
	
	
}

void WallDistanceStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ) {

	FLOAT posY = this->_parameters.meshsize->getPosY(i, j, k);
	FLOAT posZ = this->_parameters.meshsize->getPosZ(i, j, k);
	FLOAT dy2 = this->_parameters.meshsize->getDy(i, j, k) / 2.0;
	FLOAT dz2 = this->_parameters.meshsize->getDz(i, j, k) / 2.0;

	// Implementation for channel case without step only
	FLOAT minY = std::min(posY+dy2, _lengthY-posY-dy2);
	FLOAT minZ = std::min(posZ+dz2, _lengthZ-posZ-dz2);
	flowField.getWallDistance().getScalar(i, j, k) = std::abs(std::min(minY, minZ));
	
}

