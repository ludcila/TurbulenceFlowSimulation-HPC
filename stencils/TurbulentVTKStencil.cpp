#include <iostream>
#include <sstream>
#include "TurbulentVTKStencil.h"

TurbulentVTKStencil::TurbulentVTKStencil(Parameters & parameters) : VTKStencil<TurbulentFlowField>(parameters) {
}

void TurbulentVTKStencil::apply(TurbulentFlowField & flowField, int i, int j) {

	const int obstacle = flowField.getFlags().getValue(i, j);

	if(obstacle & OBSTACLE_SELF) {
		this->_pressureStringStream << "0.0" << std::endl;
		this->_velocityStringStream << "0.0 0.0 0.0" << std::endl;
		this->_turbulentViscosityStringStream << "0.0" << std::endl;
	} else {
		FLOAT pressure;
		FLOAT velocity[2];
		flowField.getPressureAndVelocity(pressure, velocity, i, j);
		this->_pressureStringStream << pressure << std::endl;
		this->_velocityStringStream << velocity[0] << " " << velocity[1] << " 0.0" << std::endl;
		this->_turbulentViscosityStringStream << flowField.getTurbulentViscosity().getScalar(i, j) << std::endl;
	}
	
}

void TurbulentVTKStencil::apply(TurbulentFlowField & flowField, int i, int j, int k) {

	const int obstacle = flowField.getFlags().getValue(i, j, k);

	if(obstacle & OBSTACLE_SELF) {
		this->_pressureStringStream << "0.0" << std::endl;
		this->_velocityStringStream << "0.0 0.0 0.0" << std::endl;
		this->_turbulentViscosityStringStream << "0.0" << std::endl;
	} else {
		FLOAT pressure;
		FLOAT velocity[3];
		flowField.getPressureAndVelocity(pressure, velocity, i, j, k);
		this->_pressureStringStream << pressure << std::endl;
		this->_velocityStringStream << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;
		this->_turbulentViscosityStringStream << flowField.getTurbulentViscosity().getScalar(i, j, k) << std::endl;
	}
	
}


void TurbulentVTKStencil::write(TurbulentFlowField & flowField, int timeStep, std::string foldername) {

	std::cout << "=== Writing VTK Output ===" << std::endl;

	// Open the file and set precision
	this->_outputFile->open(this->getFilename(timeStep, foldername).c_str());
	*this->_outputFile << std::fixed << std::setprecision(6);

	// Output the different sections of the file
	this->writeFileHeader();
	this->writeGrid(flowField);
	this->writeCellDataHeader(flowField);
	this->writePressure();
	this->writeVelocity();
	this->writeTurbulentViscosity();

	// Close the file
	this->_outputFile->close();

	// Clear string streams
	this->clearStringStreams();
	
}

void TurbulentVTKStencil::writeTurbulentViscosity() {

	// Print header
	*this->_outputFile << "SCALARS turbulent_viscosity float 1" << std::endl;
	*this->_outputFile << "LOOKUP_TABLE default" << std::endl;

	// Print pressure values
	*this->_outputFile << this->_turbulentViscosityStringStream.str().c_str();

	*this->_outputFile << std::endl;

}

void TurbulentVTKStencil::clearStringStreams() {
	_pressureStringStream.str("");
	_velocityStringStream.str("");
	_turbulentViscosityStringStream.str("");
}
