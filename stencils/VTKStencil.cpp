#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include "VTKStencil.h"
#include "Iterators.h"
#include <string>
#include <fstream>
#include <sstream>

VTKStencil::VTKStencil ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) {

}

void VTKStencil::apply ( FlowField & flowField, int i, int j ) {

}

void VTKStencil::apply ( FlowField & flowField, int i, int j, int k ) {
	if(i == 1 || j == 1 || k == 1) return;
    const int obstacle = flowField.getFlags().getValue(i, j, k);
	FLOAT pressure;
	FLOAT velocity[3];
	flowField.getPressureAndVelocity(pressure, velocity, i, j, k);
	if(obstacle & OBSTACLE_SELF) {
		*this->_outputFile << "0" << std::endl;
	} else {
		*this->_outputFile << pressure << std::endl;
	}
}

void VTKStencil::write ( FlowField & flowField, int timeStep ) {

	std::cout << "=== Writing VTK Output ===" << this->_parameters.vtk.prefix << std::endl;
	
	// Generate the filename
	std::stringstream filename;
	filename << this->_parameters.vtk.prefix << "_" << timeStep << ".vtk";
	std::cout << filename.str() << std::endl;
	
	// Open the file
	std::ofstream vtkFile;
	vtkFile.open(filename.str().c_str());
	this->_outputFile = &vtkFile;
	
	// Print file header
	vtkFile << "# vtk DataFile Version 2.0" << std::endl;
	vtkFile << "Boom" << std::endl;
	vtkFile << "ASCII" << std::endl << std::endl;
	
	// Print grid header
	int nx = flowField.getNx() + 1;
	int ny = flowField.getNy() + 1;
	int nz = flowField.getNz() + 1;
	vtkFile << "DATASET STRUCTURED_GRID" << std::endl;
	vtkFile << "DIMENSIONS " << nx << " " << ny  << " " << nz << std::endl;
	vtkFile << "POINTS " << nx*ny*nz << " float" << std::endl;
	
	// Print grid
	for(int k = 2; k < flowField.getCellsZ(); k++) {
		for(int j = 2; j < flowField.getCellsY(); j++) {
			for(int i = 2; i < flowField.getCellsX(); i++) {
				vtkFile << this->_parameters.meshsize->getPosX(i, j, k) << " " << this->_parameters.meshsize->getPosY(i, j, k) << " " << this->_parameters.meshsize->getPosZ(i, j, k) << std::endl;
			}
		}
	}
	
	// Print pressure
	vtkFile << std::endl;
	vtkFile << "CELL_DATA " << flowField.getNx()*flowField.getNy()*flowField.getNz() << std::endl;
	vtkFile << "SCALARS pressure float 1" << std::endl;
	vtkFile << "LOOKUP_TABLE default" << std::endl;
	FieldIterator<FlowField> vtkIt(flowField, this->_parameters, *this);
	vtkIt.iterate();
	
	vtkFile << std::endl;
	
	this->_outputFile->close();
	
}
