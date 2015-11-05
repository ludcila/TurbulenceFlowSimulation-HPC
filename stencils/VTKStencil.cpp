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
	if(i == 1 || j == 1) return;
    const int obstacle = flowField.getFlags().getValue(i, j);
	FLOAT pressure;
	FLOAT velocity[2];
	flowField.getPressureAndVelocity(pressure, velocity, i, j);
	if(obstacle & OBSTACLE_SELF) {
		if(this->_outputPressure) {
			*this->_outputFile << "0.0" << std::endl;
		} else {
			*this->_outputFile << "0.0 0.0 0.0" << std::endl;
		}
	} else {
		if(this->_outputPressure) {
			*this->_outputFile << pressure << std::endl;
		} else {
			*this->_outputFile << velocity[0] << " " << velocity[1] << " 0.0" << std::endl;
		}
	}
}

void VTKStencil::apply ( FlowField & flowField, int i, int j, int k ) {
	if(i == 1 || j == 1 || k == 1) return;
    const int obstacle = flowField.getFlags().getValue(i, j, k);
	FLOAT pressure;
	FLOAT velocity[3];
	flowField.getPressureAndVelocity(pressure, velocity, i, j, k);
	if(obstacle & OBSTACLE_SELF) {
		if(this->_outputPressure) {
			*this->_outputFile << "0.0" << std::endl;
		} else {
			*this->_outputFile << "0.0 0.0 0.0" << std::endl;
		}
	} else {
		if(this->_outputPressure) {
			*this->_outputFile << pressure << std::endl;
		} else {
			*this->_outputFile << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;
		}
	}
}

void VTKStencil::write ( FlowField & flowField, int timeStep ) {

	int is3D = this->_parameters.geometry.dim == 3;

	std::cout << "=== Writing VTK Output ===" << this->_parameters.vtk.prefix << std::endl;
	
	// Generate the filename
	std::stringstream filename;
	filename << this->_parameters.vtk.prefix << "." << timeStep << ".vtk";
	std::cout << filename.str() << std::endl;
	
	// Open the file
	this->_outputFile->open(filename.str().c_str());
	*this->_outputFile << std::fixed << std::setprecision(8);
	
	// Print file header
	*this->_outputFile << "# vtk DataFile Version 2.0" << std::endl;
	*this->_outputFile << "Boom" << std::endl;
	*this->_outputFile << "ASCII" << std::endl << std::endl;
	
	// Print grid header
	int nx = flowField.getNx();
	int ny = flowField.getNy();
	int nz = is3D ? flowField.getNz() : 0;
	*this->_outputFile << "DATASET STRUCTURED_GRID" << std::endl;
	*this->_outputFile << "DIMENSIONS " << (nx+1) << " " << (ny+1)  << " " << (nz+1) << std::endl;
	*this->_outputFile << "POINTS " << (nx+1)*(ny+1)*(nz+1) << " float" << std::endl;
	
	int cellsX = flowField.getCellsX();
	int cellsY = flowField.getCellsY();
	int cellsZ = is3D ? flowField.getCellsZ() : 3;
	// Print grid
	for(int k = 2; k < cellsZ; k++) {
		for(int j = 2; j < cellsY; j++) {
			for(int i = 2; i < cellsX; i++) {
				*this->_outputFile << this->_parameters.meshsize->getPosX(i, j, k) << " " << this->_parameters.meshsize->getPosY(i, j, k) << " " << this->_parameters.meshsize->getPosZ(i, j, k) << std::endl;
			}
		}
	}
	*this->_outputFile << std::endl;
	
	FieldIterator<FlowField> vtkIt(flowField, this->_parameters, *this);
	
	// Print pressure
	this->_outputPressure = 1;
	int numCells = is3D ? nx*ny*nz : nx*ny;
	*this->_outputFile << "CELL_DATA " << numCells << std::endl;
	*this->_outputFile << "SCALARS pressure float 1" << std::endl;
	*this->_outputFile << "LOOKUP_TABLE default" << std::endl;
	vtkIt.iterate();
	*this->_outputFile << std::endl;
	
	// Print velocity
	this->_outputPressure = 0;
	*this->_outputFile << "VECTORS velocity float" << std::endl;
	vtkIt.iterate();
	*this->_outputFile << std::endl;
	
	this->_outputFile->close();
	
}
