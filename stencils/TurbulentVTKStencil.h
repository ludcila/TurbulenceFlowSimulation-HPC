#include <string>
#include "VTKStencil.h"
#include "TurbulentFlowField.h"
#include "Parameters.h"

class TurbulentVTKStencil : public VTKStencil, public FieldStencil<TurbulentFlowField> {
	
	protected:
		std::stringstream _turbulentViscosityStringStream;
	public:
		TurbulentVTKStencil(Parameters & parameters);
		void apply(TurbulentFlowField & flowField, int i, int j);
		void apply(TurbulentFlowField & flowField, int i, int j, int k);
		void write(TurbulentFlowField & flowField, int timeStep, std::string foldername);
		void writeTurbulentViscosity();
		void clearStringStreams();
		
};
