#ifndef _TURBULENT_KEPS_VTK_STENCIL_H_
#define _TURBULENT_KEPS_VTK_STENCIL_H_



#include <string>
#include "TurbulentVTKStencil.h"
#include "TurbulentFlowField.h"
#include "Parameters.h"

class TurbulentKepsVTKStencil : public TurbulentVTKStencil {
	
	protected:
		std::stringstream _turbulentKineticEnergyStringStream;
		std::stringstream _turbulentDissipationRateStringStream;
	public:
		TurbulentKepsVTKStencil(Parameters & parameters);
		void apply(TurbulentFlowField & flowField, int i, int j);
		void apply(TurbulentFlowField & flowField, int i, int j, int k);
		void write(TurbulentFlowField & flowField, int timeStep, std::string foldername);
		void writeTurbulentKineticEnergy();
		void writeTurbulentDissipationRate();
		void clearStringStreams();
		
};

#endif
