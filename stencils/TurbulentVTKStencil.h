#include <string>
#include "VTKStencil.h"
#include "TurbulentFlowField.h"
#include "Parameters.h"

/*
	VTKStencil is now implemented as a template class.
	
	The compiler will generate the (default) methods specified in the h file for both
	VTKStencil<FlowField> and VTKStencil<TurbulentFlowField> (and any other type of field if specified).
	At the moment, for VTKStencil<TurbulentFlowField> we need to specialize (override) and add some methods.
	This can be done in this file.
	
	Normally, we would specialize member functions this way:
	template<>
	void VTKStencil<TurbulentFlowField>::apply(TurbulentFlowField & flowField, int i, int j, int k) {
		// code here
	}
	
	But apparently, there is no way to add a member variable without rewriting everything inside the template,
	so the only solution is to inherit from VTKStencil<TurbulentFlowField> and add the member variables and
	override the functions there :(
	(https://bytes.com/topic/c/answers/805802-adding-members-template-class-specialization)
	
*/


class TurbulentVTKStencil : public VTKStencil<TurbulentFlowField> {
	
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
