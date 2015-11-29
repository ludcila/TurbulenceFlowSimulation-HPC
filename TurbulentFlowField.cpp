#include "TurbulentFlowField.h"

TurbulentFlowField::TurbulentFlowField ( int Nx, int Ny ) :
	FlowField(Nx, Ny),
	_nitau( ScalarField ( Nx + 3, Ny + 3 ) ),
	_distance( ScalarField ( Nx + 3, Ny + 3 ) ) {}


TurbulentFlowField::TurbulentFlowField ( int Nx, int Ny, int Nz ) :
	FlowField(Nx, Ny, Nz),
	_nitau( ScalarField ( Nx + 3, Ny + 3, Nz+3 ) ),
	_distance( ScalarField ( Nx + 3, Ny + 3, Nz+3 ) ) {}


TurbulentFlowField::TurbulentFlowField (const Parameters & parameters):
	FlowField(parameters),
	_nitau(parameters.geometry.dim==2?ScalarField(_size_x + 3, _size_y + 3):
                      ScalarField(_size_x + 3, _size_y + 3, _size_z + 3)),
	_distance(parameters.geometry.dim==2?ScalarField(_size_x + 3, _size_y + 3):
                      ScalarField(_size_x + 3, _size_y + 3, _size_z + 3)) {}

ScalarField & TurbulentFlowField::getTurbulentViscosity () {
    return _nitau;
}