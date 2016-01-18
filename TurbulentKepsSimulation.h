#ifndef _TURBULENT_KEPS_SIMULATION_H_
#define _TURBULENT_KEPS_SIMULATION_H_

#include "Simulation.h"
#include "TurbulentFlowField.h"
#include "stencils/TurbulentViscosityKepsStencil.h"
#include "stencils/TurbulentKepsVTKStencil.h"
#include "stencils/TurbulenceFGHStencil.h"
#include "stencils/TurbulentKepsStencil.h"
#include "stencils/KepsBoundaryStencil.h"
#include "stencils/MaxNuStencil.h"
#include "stencils/FmuStencil.h"
#include "parallelManagers/PetscParallelManagerTurbulent.h"

class TurbulentKepsSimulation : public Simulation {

	protected:
		TurbulentFlowField &_turbulentFlowField;
    	FieldIterator<TurbulentFlowField> _turbulentFghIterator;
		TurbulenceFGHStencil _turbulentFghStencil; // K-eps
    	FieldIterator<TurbulentFlowField> _turbulentVtkIterator;
		TurbulentKepsVTKStencil _turbulentVtkStencil;
		TurbulentViscosityKepsStencil _turbulentViscosityStencil; // Viscosity in K-eps
		TurbulentKepsStencil _turbulentKepsStencil; // K and eps
		FieldIterator<TurbulentFlowField> _turbulentViscosityIterator;
		FieldIterator<TurbulentFlowField> _turbulentKepsIterator;
		MaxNuStencil _maxNuStencil;
		FieldIterator<TurbulentFlowField> _maxNuFieldIterator;
        GlobalBoundaryIterator<TurbulentFlowField> _maxNuBoundaryIterator;
        KepsBoundaryStencil _kepsBoundaryStencil;
        GlobalBoundaryIterator<TurbulentFlowField> _kepsBoundaryIterator;
        //GlobalBoundaryIterator<TurbulentFlowField> _turbulentViscosityBoundaryIterator;
		//TurbulentViscosityBoundaryStencil _turbulentViscosityBoundaryStencil;
		FieldIterator<TurbulentFlowField> _fmuIterator;
		FmuStencil _fmuStencil;
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
