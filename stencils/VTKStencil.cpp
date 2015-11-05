#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include "VTKStencil.h"
#include "Iterators.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

VTKStencil::VTKStencil ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) {
	this->_outputFile = new std::ofstream;
}

VTKStencil::~VTKStencil () {
	delete this->_outputFile;
}

void VTKStencil::apply ( FlowField & flowField, int i, int j ) {
	
	// Skip the ghost cell
	if(i == 1 || j == 1) return;
	
    const int obstacle = flowField.getFlags().getValue(i, j);
    
	if(obstacle & OBSTACLE_SELF) {
		if(this->isWritingPressure()) {
			*this->_outputFile << "0.0" << std::endl;
		} else {
			*this->_outputFile << "0.0 0.0 0.0" << std::endl;
		}
	} else {
		FLOAT pressure;
		FLOAT velocity[2];
		flowField.getPressureAndVelocity(pressure, velocity, i, j);
		if(this->isWritingPressure()) {
			*this->_outputFile << pressure << std::endl;
		} else {
			*this->_outputFile << velocity[0] << " " << velocity[1] << " 0.0" << std::endl;
		}
	}
	
}

void VTKStencil::apply ( FlowField & flowField, int i, int j, int k ) {

	// Skip the ghost cell
	if(i == 1 || j == 1 || k == 1) return;
	
    const int obstacle = flowField.getFlags().getValue(i, j, k);
    
	if(obstacle & OBSTACLE_SELF) {
		if(this->isWritingPressure()) {
			*this->_outputFile << "0.0" << std::endl;
		} else {
			*this->_outputFile << "0.0 0.0 0.0" << std::endl;
		}
	} else {
		FLOAT pressure;
		FLOAT velocity[3];
		flowField.getPressureAndVelocity(pressure, velocity, i, j, k);
		if(this->isWritingPressure()) {
			*this->_outputFile << pressure << std::endl;
		} else {
			*this->_outputFile << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;
		}
	}
	
}

void VTKStencil::write ( FlowField & flowField, int timeStep, std::string foldername ) {

	std::cout << "=== Writing VTK Output ===" << std::endl;
	
	// Open the file and set precision
	this->_outputFile->open(this->getFilename(timeStep, foldername).c_str());
	*this->_outputFile << std::fixed << std::setprecision(6);
	
	// Output the different sections of the file
	this->writeFileHeader();
	this->writeGrid(flowField);
	this->writeCellDataHeader(flowField);
	this->writePressure(flowField);
	this->writeVelocity(flowField);
	
	// Close the file
	this->_outputFile->close();
	
}

void VTKStencil::writeFileHeader() {
	*this->_outputFile << "# vtk DataFile Version 2.0" << std::endl;
	*this->_outputFile << "Turbulent Flow Simulation" << std::endl;
	*this->_outputFile << "ASCII" << std::endl << std::endl;
}

void VTKStencil::writeGrid ( FlowField & flowField ) {

	int is3D = this->_parameters.geometry.dim == 3;
	
	// Get number of grid points
	int gridX = flowField.getNx() + 1;
	int gridY = flowField.getNy() + 1;
	int gridZ = is3D ? flowField.getNz() + 1 : 1;
	
	// Print grid header
	*this->_outputFile << "DATASET STRUCTURED_GRID" << std::endl;
	*this->_outputFile << "DIMENSIONS " << gridX << " " << gridY  << " " << gridZ << std::endl;
	*this->_outputFile << "POINTS " << gridX * gridY * gridZ << " float" << std::endl;
	
	// Get id of first and last cell in each dimension
	int firstCell = 2; // 0 and 1 are ghost cells
	int lastCellX = flowField.getCellsX();
	int lastCellY = flowField.getCellsY();
	int lastCellZ = is3D ? flowField.getCellsZ() : 3;
	
	// Print grid
	for(int k = firstCell; k < lastCellZ; k++) {
		for(int j = firstCell; j < lastCellY; j++) {
			for(int i = firstCell; i < lastCellX; i++) {
				*this->_outputFile << this->_parameters.meshsize->getPosX(i, j, k) << " " << this->_parameters.meshsize->getPosY(i, j, k) << " " << this->_parameters.meshsize->getPosZ(i, j, k) << std::endl;
			}
		}
	}
	
	*this->_outputFile << std::endl;
	
}

void VTKStencil::writeCellDataHeader ( FlowField & flowField ) {

	int is3D = this->_parameters.geometry.dim == 3;
	int cellsX = flowField.getNx();
	int cellsY = flowField.getNy();
	int cellsZ = flowField.getNz();
	int numCells = is3D ? (cellsX * cellsY * cellsZ) : (cellsX * cellsY);
	*this->_outputFile << "CELL_DATA " << numCells << std::endl;
	
}

void VTKStencil::writePressure ( FlowField & flowField ) {

	// Print header
	*this->_outputFile << "SCALARS pressure float 1" << std::endl;
	*this->_outputFile << "LOOKUP_TABLE default" << std::endl;
	
	// Set flag to output pressure
	this->isWritingPressure(1);
	
	// Iterate over the flow field
	FieldIterator<FlowField> vtkIt(flowField, this->_parameters, *this);
	vtkIt.iterate();
	
	*this->_outputFile << std::endl;
	
}

void VTKStencil::writeVelocity ( FlowField & flowField ) {

	// Print header
	*this->_outputFile << "VECTORS velocity float" << std::endl;

	// Set flag to output velocity
	this->isWritingVelocity(1);
	
	// Iterate over the flow field
	FieldIterator<FlowField> vtkIt(flowField, this->_parameters, *this);
	vtkIt.iterate();
	
	*this->_outputFile << std::endl;
	
}

std::string VTKStencil::getFilename( int timeStep, std::string foldername ) {
	std::stringstream filename;
	filename << foldername << "/" << this->_parameters.vtk.prefix << "." << timeStep << ".vtk";
	return filename.str();
}

void VTKStencil::isWritingPressure(bool isWritingPressure) {
	this->_isWritingPressure = isWritingPressure;
}
void VTKStencil::isWritingVelocity(bool isWritingVelocity) {
	this->_isWritingPressure = !isWritingVelocity;
}
bool VTKStencil::isWritingPressure() {
	return this->_isWritingPressure;
}
bool VTKStencil::isWritingVelocity() {
	return !this->_isWritingPressure;
}


