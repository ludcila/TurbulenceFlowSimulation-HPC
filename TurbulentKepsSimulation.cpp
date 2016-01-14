#include "TurbulentKepsSimulation.h"
#include "stencils/WallDistanceStencil.h"

TurbulentKepsSimulation::TurbulentKepsSimulation(Parameters &parameters, TurbulentFlowField &flowField):
	Simulation(parameters, flowField),
	_turbulentFlowField(flowField),
	_turbulentFghStencil(parameters),
	_turbulentFghIterator(_turbulentFlowField, parameters, _turbulentFghStencil),
	_turbulentViscosityStencil(parameters),
	_turbulentViscosityIterator(_turbulentFlowField, parameters, _turbulentViscosityStencil, 1, 0),
	_turbulentKepsStencil(parameters),
	_turbulentKepsIterator(_turbulentFlowField, parameters, _turbulentKepsStencil, 1, 0),
	_turbulentVtkStencil(parameters),
	_turbulentVtkIterator(_turbulentFlowField, parameters, _turbulentVtkStencil, 1, 0),
	_kepsBoundaryStencil(parameters),
	_kepsBoundaryIterator(_turbulentFlowField, parameters, _kepsBoundaryStencil, 1, 0),
	_maxNuStencil(parameters),
    _maxNuFieldIterator(_turbulentFlowField,parameters,_maxNuStencil),
    _maxNuBoundaryIterator(_turbulentFlowField,parameters,_maxNuStencil),
    _fmuStencil(parameters),
    _fmuIterator(_turbulentFlowField, parameters,_fmuStencil,1,0),
	_parallelManagerTurbulent(_turbulentFlowField, parameters)
{
}

void TurbulentKepsSimulation::solveTimestep(){

	// determine and set max. timestep which is allowed in this simulation
	setTimeStep();
	
	// compute k, eps, and viscosity
std::cout<<"------------------------------------"<<std::endl;
	_turbulentFlowField.getFmu().show();
	_fmuIterator.iterate();
	_kepsBoundaryIterator.iterate();
	_turbulentKepsIterator.iterate();
	_turbulentViscosityIterator.iterate();
	_parallelManagerTurbulent.communicateViscosity();
	
	// compute fgh
	_turbulentFghIterator.iterate();
	_wallFGHIterator.iterate();
	
	// compute the right hand side
	_rhsIterator.iterate();
	
	// solve for pressure
	_solver.solve();
	
	// communicate pressure values
	_parallelManagerTurbulent.communicatePressure();
	
	// compute velocity
	_velocityIterator.iterate();
	
	// set obstacle boundaries
	_obstacleIterator.iterate();
	
	// communicate velocity values
	_parallelManagerTurbulent.communicateVelocity();
	
	// Iterate for velocities on the boundary
	_wallVelocityIterator.iterate();
	
}

void TurbulentKepsSimulation::setTimeStep(){

	const FLOAT cinematicviscosity=1.0/_parameters.flow.Re;
	FLOAT localMin, globalMin;
	assertion(_parameters.geometry.dim == 2 || _parameters.geometry.dim == 3);
	FLOAT factor = 1.0/(_parameters.meshsize->getDxMin() * _parameters.meshsize->getDxMin()) +
		         1.0/(_parameters.meshsize->getDyMin() * _parameters.meshsize->getDyMin());

	// determine maximum velocity
	_maxUStencil.reset();
	_maxUFieldIterator.iterate();
	_maxUBoundaryIterator.iterate();

	_maxNuStencil.reset();
	_maxNuFieldIterator.iterate();
	_maxNuBoundaryIterator.iterate();

	if (_parameters.geometry.dim == 3) {
	factor += 1.0/(_parameters.meshsize->getDzMin() * _parameters.meshsize->getDzMin());
	_parameters.timestep.dt = 1.0 / _maxUStencil.getMaxValues()[2];
	} else {
	_parameters.timestep.dt = 1.0 / _maxUStencil.getMaxValues()[0];
	}

	localMin = std::min(_parameters.timestep.dt,
		                            std::min(std::min(1.0/_maxNuStencil.getMaxValue(),
		                            1.0 / _maxUStencil.getMaxValues()[0]),
		                            1.0 / _maxUStencil.getMaxValues()[1]));

	// Here, we select the type of operation before compiling. This allows to use the correct
	// data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
	// machines.

	globalMin = MY_FLOAT_MAX;
	MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

	_parameters.timestep.dt = globalMin;
	_parameters.timestep.dt *= _parameters.timestep.tau;

}

void TurbulentKepsSimulation::plotVTK(int timeStep, std::string foldername) {
	_turbulentVtkIterator.iterate();
	_turbulentVtkStencil.write(this->_turbulentFlowField, timeStep, foldername);
}

void TurbulentKepsSimulation::initializeFlowField() {
	Simulation::initializeFlowField();
	WallDistanceStencil wds(_parameters);
	FieldIterator<TurbulentFlowField> it(_turbulentFlowField, _parameters, wds);
	it.iterate();
	FLOAT eps0 = _turbulentFlowField.getDissipationRate().getScalar(1, 2);
	eps0 = 1; // hardcode initialization for now
    if (_parameters.geometry.dim==2){
		const int sizex = _flowField.getNx();
		const int sizey = _flowField.getNy();
		for (int i =1 ;i < sizex+3;i++) {
			for (int j =1 ;i < sizey+3;i++) {
				_turbulentFlowField.getDissipationRate().getScalar(i,j) = eps0;
			}
		}
    } else {
		const int sizex = _flowField.getNx();
		const int sizez = _flowField.getNz();
		const int sizey = _flowField.getNy();
    }
};
