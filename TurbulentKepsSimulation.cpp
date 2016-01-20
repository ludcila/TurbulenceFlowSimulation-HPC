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
	_kepsObstacleStencil(parameters),
	_kepsObstacleIterator(_turbulentFlowField, parameters, _kepsObstacleStencil),
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
	_fmuIterator.iterate();
	_kepsObstacleIterator.iterate();
	_kepsBoundaryIterator.iterate();
	_turbulentKepsIterator.iterate();
	_turbulentFlowField.swapKeps();
	_kepsBoundaryIterator.iterate();
	
	_parallelManagerTurbulent.communicateKineticEnergy();
	_parallelManagerTurbulent.communicateDissipationRate();
	
	//_fmuIterator.iterate();
	_turbulentViscosityIterator.iterate();
	_kepsBoundaryIterator.iterate();
	
	
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
	
	// communicate velocity values
	_parallelManagerTurbulent.communicateVelocity();
	
	// Iterate for velocities on the boundary
	_wallVelocityIterator.iterate();
	
	// set obstacle boundaries
	_obstacleIterator.iterate();

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
	
	// To be moved to somewhere else, in case the formula works
	// _parameters.turbulence.gamma = _parameters.timestep.dt * std::max(_maxUStencil.getMaxValues()[0], _maxUStencil.getMaxValues()[1]);
	// _parameters.turbulence.gamma = 1;

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
	
	FLOAT kin =1.5* pow(_parameters.walls.vectorLeft[0]*0.05,2);
	FLOAT epsin;

	if (_parameters.geometry.dim==2)
	 epsin= 0.09/(0.038*_parameters.geometry.lengthY)*pow(kin,1.5);

	else
	 epsin= 0.09/(0.038*_parameters.geometry.lengthY*_parameters.geometry.lengthZ*2/(_parameters.geometry.lengthY+_parameters.geometry.lengthZ))*pow(kin,1.5);
	
    if (_parameters.geometry.dim==2){
		const int sizex = _flowField.getNx();
		const int sizey = _flowField.getNy();
		for (int i =1 ;i < sizex+3; i++) {
			for (int j =1 ;j < sizey+3; j++) {
				_turbulentFlowField.getDissipationRate().getScalar(i, j) = epsin;
				_turbulentFlowField.getKineticEnergy().getScalar(i, j) = kin;
			}
		}
    } else {
		const int sizex = _flowField.getNx();
		const int sizey = _flowField.getNy();
		const int sizez = _flowField.getNz();
	
		for (int i =1 ;i < sizex+3; i++) {
			for (int j =1 ;j < sizey+3; j++) {
				for(int k=1; k<sizez+3; k++){
					_turbulentFlowField.getDissipationRate().getScalar(i, j, k) = epsin;
					_turbulentFlowField.getKineticEnergy().getScalar(i, j, k) = kin;
				}
			}
		}
    }
    
    // Initialize turbulent viscosity by iterating it once
	_fmuIterator.iterate();
	_turbulentViscosityIterator.iterate();
	
};
