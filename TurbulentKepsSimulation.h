#ifndef _TURBULENT_KEPS_SIMULATION_H_
#define _TURBULENT_KEPS_SIMULATION_H_

#include "Simulation.h"
#include "TurbulentFlowField.h"
#include "stencils/TurbulentViscosityKepsStencil.h"
#include "stencils/TurbulentVTKStencil.h"
#include "stencils/TurbulenceFGHStencil.h"
#include "stencils/MaxNuStencil.h"
#include "parallelManagers/PetscParallelManagerTurbulent.h"

class TurbulentKepsSimulation : public Simulation {

	protected:
		TurbulentFlowField &_turbulentFlowField;
    	FieldIterator<TurbulentFlowField> _turbulentFghIterator;
		TurbulenceFGHStencil _turbulentFghStencil; // K-eps
    	FieldIterator<TurbulentFlowField> _turbulentVtkIterator;
		TurbulentVTKStencil _turbulentVtkStencil;
		TurbulentViscosityKepsStencil _turbulentViscosityStencil; // K-eps
		FieldIterator<TurbulentFlowField> _turbulentViscosityIterator;
		MaxNuStencil _maxNuStencil;
		FieldIterator<TurbulentFlowField> _maxNuFieldIterator;
        GlobalBoundaryIterator<TurbulentFlowField> _maxNuBoundaryIterator;
        //GlobalBoundaryIterator<TurbulentFlowField> _turbulentViscosityBoundaryIterator;
		//TurbulentViscosityBoundaryStencil _turbulentViscosityBoundaryStencil;
		PetscParallelManagerTurbulent _parallelManagerTurbulent;

	public:

		TurbulentKepsSimulation(Parameters &parameters, TurbulentFlowField &flowField);
		virtual ~TurbulentKepsSimulation(){}
		virtual void solveTimestep();
		virtual void plotVTK(int timeStep, std::string foldername);
		virtual void initializeFlowField();


	protected:
		virtual void setTimeStep();

};

#endif // _TURBULENT_SIMULATION_H_
