#include <iostream>
#include <sstream>
#include <iomanip>
#include "TurbulentKepsVTKStencil.h"
#include "StencilFunctions.h"

TurbulentKepsVTKStencil::TurbulentKepsVTKStencil(Parameters & parameters) :
	TurbulentVTKStencil(parameters)
{
}

void TurbulentKepsVTKStencil::apply(TurbulentFlowField & flowField, int i, int j) {

	TurbulentVTKStencil::apply(flowField, i, j);

	const int obstacle = flowField.getFlags().getValue(i, j);

	if(obstacle & OBSTACLE_SELF) {
		_turbulentKineticEnergyStringStream << "0.0" << std::endl;
		_turbulentDissipationRateStringStream << "0.0" << std::endl;
	} else {
		_turbulentKineticEnergyStringStream << flowField.getKineticEnergy().getScalar(i, j) << std::endl;
		_turbulentDissipationRateStringStream << flowField.getDissipationRate().getScalar(i, j) << std::endl;
	}
	
}

void TurbulentKepsVTKStencil::apply(TurbulentFlowField & flowField, int i, int j, int k) {

	
}


void TurbulentKepsVTKStencil::write(TurbulentFlowField & flowField, int timeStep, std::string foldername) {

	std::cout << "=== Writing VTK Output ===" << std::endl;

	// Open the file and set precision
	_outputFile->open(getFilename(timeStep, foldername).c_str());
	*_outputFile << std::fixed << std::setprecision(6);

	// Output the different sections of the file
	writeFileHeader();
	writeGrid(flowField);
	writeCellDataHeader(flowField);
	writePressure();
	writeVelocity();
	writeTurbulentViscosity();
	writeShearStresses();
	writeTurbulentKineticEnergy();
	writeTurbulentDissipationRate();

	// Close the file
	_outputFile->close();

	// Clear string streams
	clearStringStreams();
	
}

void TurbulentKepsVTKStencil::writeTurbulentKineticEnergy() {

	// Print header
	*_outputFile << "SCALARS turbulent_k float 1" << std::endl;
	*_outputFile << "LOOKUP_TABLE default" << std::endl;

	// Print turbulent k
	*_outputFile << _turbulentKineticEnergyStringStream.str().c_str();

	*_outputFile << std::endl;

}


void TurbulentKepsVTKStencil::writeTurbulentDissipationRate() {

	// Print header
	*_outputFile << "SCALARS turbulent_eps float 1" << std::endl;
	*_outputFile << "LOOKUP_TABLE default" << std::endl;

	// Print turbulent eps
	*_outputFile << _turbulentDissipationRateStringStream.str().c_str();

	*_outputFile << std::endl;

}


void TurbulentKepsVTKStencil::clearStringStreams() {
	_pressureStringStream.str("");
	_velocityStringStream.str("");
	_turbulentViscosityStringStream.str("");
	_turbulentKineticEnergyStringStream.str("");
	_turbulentDissipationRateStringStream.str("");
	_viscousStressXYStringStream.str("");
	_viscousStressXZStringStream.str("");
	_ReStressXYStringStream.str("");
	_ReStressXZStringStream.str("");
}


