#ifndef _TURBULENT_SIMULATION_H_
#define _TURBULENT_SIMULATION_H_

#include "Simulation.h"
#include "TurbulentFlowField.h"


class TurbulentSimulation : public Simulation {
	
	protected:
		TurbulentFlowField &_flowField;

	public:
		TurbulentSimulation(Parameters &parameters, TurbulentFlowField &flowField);
		virtual ~TurbulentSimulation(){}
		virtual void solveTimestep();

	protected:
		virtual void setTimeStep();
    
};

#endif // _TURBULENT_SIMULATION_H_

