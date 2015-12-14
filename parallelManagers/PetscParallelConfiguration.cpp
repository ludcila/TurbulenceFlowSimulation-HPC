#include "PetscParallelConfiguration.h"

PetscParallelConfiguration::PetscParallelConfiguration(Parameters & parameters):
    _parameters(parameters) {

    // Obtain the rank of the current processor
    int rank, nproc;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);

    parameters.parallel.rank = rank;

    // Obtain the position of this subdomain, and locate its neighbors.
    createIndices();
    locateNeighbors();
    computeSizes();

    int nprocFromFile = _parameters.parallel.numProcessors[0] *
                        _parameters.parallel.numProcessors[1];

    if (parameters.geometry.dim == 3) {
        nprocFromFile *= _parameters.parallel.numProcessors[2];
    }

    if (nproc != nprocFromFile){
        handleError(1, "The number of processors specified in the configuration file doesn't match the communicator");
    }


    MPI_Comm_split(PETSC_COMM_WORLD, _parameters.parallel.firstCorner[0], rank, &_parameters.parallel.planeComm);
    MPI_Comm_rank(_parameters.parallel.planeComm, &parameters.parallel.plane_rank);

    set_centerline_flags();
    centerline();


}


PetscParallelConfiguration::~PetscParallelConfiguration(){
    freeSizes();
}


void PetscParallelConfiguration::locateNeighbors(){

    int i = _parameters.parallel.indices[0];
    int j = _parameters.parallel.indices[1];
    int k = _parameters.parallel.indices[2];

    if (_parameters.geometry.dim == 2){
        _parameters.parallel.leftNb   = computeRankFromIndices(_parameters, i-1, j, 0);
        _parameters.parallel.rightNb  = computeRankFromIndices(_parameters, i+1, j, 0);
        _parameters.parallel.bottomNb = computeRankFromIndices(_parameters, i, j-1, 0);
        _parameters.parallel.topNb    = computeRankFromIndices(_parameters, i, j+1, 0);

        // The following two are not used in this case
        _parameters.parallel.frontNb = MPI_PROC_NULL;
        _parameters.parallel.backNb  = MPI_PROC_NULL;
    } else {
        _parameters.parallel.leftNb   = computeRankFromIndices(_parameters, i-1, j, k);
        _parameters.parallel.rightNb  = computeRankFromIndices(_parameters, i+1, j, k);
        _parameters.parallel.bottomNb = computeRankFromIndices(_parameters, i, j-1, k);
        _parameters.parallel.topNb    = computeRankFromIndices(_parameters, i, j+1, k);
        _parameters.parallel.frontNb  = computeRankFromIndices(_parameters, i, j, k-1);
        _parameters.parallel.backNb   = computeRankFromIndices(_parameters, i, j, k+1);
    }

    // If periodic boundaries declared, let the process itself deal with it, without communication
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
}


void PetscParallelConfiguration::createIndices(){
    int & rank = _parameters.parallel.rank;
    _parameters.parallel.indices[0] = rank % _parameters.parallel.numProcessors[0];
    _parameters.parallel.indices[1] = (rank / _parameters.parallel.numProcessors[0]) %
                                      _parameters.parallel.numProcessors[1];
    _parameters.parallel.indices[2] = rank / (_parameters.parallel.numProcessors[0] *
                                              _parameters.parallel.numProcessors[1]);
}


int PetscParallelConfiguration::computeRankFromIndices(Parameters & parameters, int i, int j, int k) {

    if( i < 0 || i >= parameters.parallel.numProcessors[0] ||
        j < 0 || j >= parameters.parallel.numProcessors[1] ||
        k < 0 || k >= parameters.parallel.numProcessors[2] ){
        return MPI_PROC_NULL;
    }
    int nrank = i + j*parameters.parallel.numProcessors[0];
    if (parameters.geometry.dim == 3){
        nrank += k * parameters.parallel.numProcessors[0] * parameters.parallel.numProcessors[1];
    }
    return nrank;
}


void PetscParallelConfiguration::computeSizes(){

    int dim = _parameters.geometry.dim;

    for (int i = 0; i < dim; i++){
        _parameters.parallel.sizes[i] = new PetscInt[_parameters.parallel.numProcessors[i]];
    }

    int geometrySizes[3];
    geometrySizes[0] = _parameters.geometry.sizeX;
    geometrySizes[1] = _parameters.geometry.sizeY;
    geometrySizes[2] = _parameters.geometry.sizeZ;

    for (int i = 0; i < dim; i++){
        for (int j = 0; j < _parameters.parallel.numProcessors[i]; j++){
            _parameters.parallel.sizes[i][j] =
                geometrySizes[i] / _parameters.parallel.numProcessors[i];
            if (j < geometrySizes[i] % _parameters.parallel.numProcessors[i]){
                _parameters.parallel.sizes[i][j] ++;
            }
        }
    }

    // Locate the position of the first element of the subdomain. Useful for plotting later on.
    for (int i = 0; i < dim; i++){
        _parameters.parallel.firstCorner[i] = 0;
        for (int j = 0; j < _parameters.parallel.indices[i]; j++){
            _parameters.parallel.firstCorner[i] += _parameters.parallel.sizes[i][j];
        }
    }
    if (dim == 2){
        _parameters.parallel.firstCorner[2] = 0;
    }

    // Select the local sizes from the already computed sizes
    for (int i = 0; i < dim; i++){
        _parameters.parallel.localSize[i] =
            _parameters.parallel.sizes[i][_parameters.parallel.indices[i]];
    }

    // If the domain lies on an edge, add one to that direction, for the artificial external
    // pressures in the PETSc solver
    for (int i = 0; i < dim; i++){
        _parameters.parallel.sizes[i][0] ++;
        _parameters.parallel.sizes[i][_parameters.parallel.numProcessors[i]-1] ++;
    }
}


void PetscParallelConfiguration::freeSizes(){

    int dim = _parameters.geometry.dim;

    for (int i = 0; i < dim; i++){
        delete[] _parameters.parallel.sizes[i];
    }
}

void PetscParallelConfiguration::set_centerline_flags()
{
  int center_processor [3]={};

  _parameters.parallel.centerlineFlag=0;

  _parameters.parallel.local_center_line_index[0]=0;
  _parameters.parallel.local_center_line_index[1]=0;
  _parameters.parallel.local_center_line_index[2]=0;

  //determine the center processor in X direction
  center_processor[0]=_parameters.parallel.indices[0];

  //determine the center processor in Y direction
  if (_parameters.parallel.indices[0] < _parameters.bfStep.xRatio * _parameters.parallel.numProcessors[0] ) {
    center_processor[1] = floor ( ( _parameters.parallel.numProcessors[1] * ( 1 + _parameters.bfStep.yRatio ) ) / 2 );
  }
  else
  {
    center_processor[1] = floor ( ( _parameters.geometry.lengthY / 2 ) / ( _parameters.geometry.lengthY / _parameters.parallel.numProcessors[1] ) );
  }
  _parameters.parallel.centerProcessor[1]=center_processor[1];

  //determine the center processor in Z direction
  center_processor[2] = floor ( _parameters.parallel.numProcessors[2] / 2.0 );
  _parameters.parallel.centerProcessor[2]=center_processor[2];

  //determie the local index of the centerline in the center processor
  if ( ( center_processor[1]==_parameters.parallel.indices[1] ) && ( center_processor[2]==_parameters.parallel.indices[2] ) )
  {
    _parameters.parallel.centerlineFlag=1;

    //in Y direction
    if ( ( (_parameters.parallel.numProcessors[1] + lrint( _parameters.bfStep.yRatio * _parameters.parallel.numProcessors[1] ) ) % 2 ) == 1 )
    {
      _parameters.parallel.local_center_line_index[1]=_parameters.parallel.localSize[1]/2;
    }
    else
    {
      _parameters.parallel.local_center_line_index[1]=0;
    }

    //in Z direction
    if ( ( _parameters.parallel.numProcessors[2] % 2 ) == 1 )
    {
      _parameters.parallel.local_center_line_index[2]=_parameters.parallel.localSize[2]/2;
    }
    else
    {
      _parameters.parallel.local_center_line_index[2]=0;
    }

    //std::cout << std::endl << "##########################" <<  _parameters.parallel.localSize[0] << "#########################" << std::endl;
    //std::cout << std::endl << "##########################" <<  _parameters.parallel.firstCorner[0] << "#########################" << std::endl;
    //std::cout << std::endl << "############oodle.tum.de/##############" <<  this->_parameters.meshsize->getPosY(2, 2) << "#########################" << std::endl;



  }

}


void PetscParallelConfiguration::centerline(){

  _parameters.parallel.centerlineFlag=0;

  int * centerLineX = new int[_parameters.parallel.numProcessors[0]];
  int * centerLineY = new int[_parameters.parallel.numProcessors[0]];
  int * centerLineZ = new int[_parameters.parallel.numProcessors[0]];


  //  x direction
  centerLineX[0]=0;
  //std::cout << std::endl << "########################## " << centerLineX[0] << " #########################" << std::endl;
  //std::cout << std::endl << "########################## " << _parameters.parallel.firstCorner[0] << " #########################" << std::endl;
  for (size_t i = 1; i < _parameters.parallel.numProcessors[0]; i++) {
    centerLineX[i]=centerLineX[i-1]+_parameters.parallel.localSize[0];
    //std::cout << std::endl << "########################## " << _parameters.parallel.firstCorner[i] << " #########################" << std::endl;
    //std::cout << std::endl << "########################## " << centerLineX[i] << " #########################" << std::endl;
  }

  //  y direction
  for (size_t i = 0; i < _parameters.parallel.numProcessors[0]; i++) {

    if (_parameters.parallel.numProcessors[1]==1)
      {
        centerLineY[i] = _parameters.geometry.sizeY/2;
      }
    else
      {
        if (_parameters.parallel.firstCorner[0] < _parameters.bfStep.xRatio * _parameters.geometry.sizeX )
          {
            centerLineY[i] = (1 + _parameters.bfStep.yRatio) * _parameters.geometry.sizeY / 2;
          }
        else
          {
            centerLineY[i] = _parameters.geometry.sizeY / 2;
          }
      }
    //std::cout << std::endl << "########################## " << centerLineY[i] << " #########################" << std::endl;
  }

    //  z direction
  for (size_t i = 0; i < _parameters.parallel.numProcessors[0]; i++) {
    centerLineZ[i] = _parameters.geometry.sizeZ / 2;
    //std::cout << std::endl << "########################## " << centerLineZ[i] << " #########################" << std::endl;
    //std::cout << std::endl << "########################## " << _parameters.parallel.firstCorner[2] << " ## " << centerLineZ[i] << " ##################" << std::endl;
  }


  if (_parameters.geometry.dim == 2) {
    for (size_t i = 0; i < _parameters.parallel.numProcessors[0]; i++) {

      if ( centerLineX[i] == _parameters.parallel.firstCorner[0] &&
           centerLineY[i] == _parameters.parallel.firstCorner[1] )
           //centerLineZ[i] == _parameters.parallel.firstCorner[2]   )
      {
        _parameters.parallel.centerlineFlag=1;
        _parameters.parallel.plane_root=_parameters.parallel.plane_rank;

        if ( ( (_parameters.parallel.numProcessors[1] + lrint( _parameters.bfStep.yRatio * _parameters.parallel.numProcessors[1] ) ) % 2 ) == 1 )
        {
          _parameters.parallel.local_center_line_index[1]=_parameters.parallel.localSize[1]/2;
        }
        else
        {
          _parameters.parallel.local_center_line_index[1]=0;
        }
          //in Z direction
        if ( ( _parameters.parallel.numProcessors[2] % 2 ) == 1 )
        {
          _parameters.parallel.local_center_line_index[2]=_parameters.parallel.localSize[2]/2;
        }
        else
        {
          _parameters.parallel.local_center_line_index[2]=0;
        }

        for (size_t i = 0; i < _parameters.parallel.numProcessors[1]; i++) {
          if (i!=_parameters.parallel.plane_rank) {

            //std::cout << std::endl << "########################## " << "inside the send loop" << " #########################" << std::endl;
            MPI_Send(&_parameters.parallel.plane_rank, 1, MPI_INT , i, 0,
                      _parameters.parallel.planeComm);
          }
        }
      }
    }

    if (!_parameters.parallel.centerlineFlag) {
      //std::cout << std::endl << "########################## " << "inside the recieve " << " #########################" << std::endl;
      MPI_Recv(&_parameters.parallel.plane_root, 1, MPI_INT , MPI_ANY_SOURCE, 0,
                _parameters.parallel.planeComm,  MPI_STATUS_IGNORE);
    }

  }
  else
  {
    for (size_t i = 0; i < _parameters.parallel.numProcessors[0]; i++) {

      if ( centerLineX[i] == _parameters.parallel.firstCorner[0] &&
           centerLineY[i] == _parameters.parallel.firstCorner[1] &&
           centerLineZ[i] == _parameters.parallel.firstCorner[2]   )
      {

        std::cout << std::endl << "########################## " << _parameters.parallel.firstCorner[1] << " ## " << centerLineY[i] << " ##################" << std::endl;

        _parameters.parallel.centerlineFlag=1;
        _parameters.parallel.plane_root=_parameters.parallel.plane_rank;

        if ( ( (_parameters.parallel.numProcessors[1] + lrint( _parameters.bfStep.yRatio * _parameters.parallel.numProcessors[1] ) ) % 2 ) == 1 )
        {
          _parameters.parallel.local_center_line_index[1]=_parameters.parallel.localSize[1]/2;
        }
        else
        {
          _parameters.parallel.local_center_line_index[1]=0;
        }
          //in Z direction
        if ( ( _parameters.parallel.numProcessors[2] % 2 ) == 1 )
        {
          _parameters.parallel.local_center_line_index[2]=_parameters.parallel.localSize[2]/2;
        }
        else
        {
          _parameters.parallel.local_center_line_index[2]=0;
        }

        for (size_t i = 0; i < _parameters.parallel.numProcessors[1]*_parameters.parallel.numProcessors[2]; i++) {
          if (i!=_parameters.parallel.plane_rank) {

//            std::cout << std::endl << "########################## " << "inside the send loop" << " #########################" << std::endl;
            MPI_Send(&_parameters.parallel.plane_rank, 1, MPI_INT , i, 0,
                      _parameters.parallel.planeComm);
          }
        }
      }
    }

    if (!_parameters.parallel.centerlineFlag) {
//      std::cout << std::endl << "########################## " << "inside the recieve " << " #########################" << std::endl;
      MPI_Recv(&_parameters.parallel.plane_root, 1, MPI_INT , MPI_ANY_SOURCE, 0,
                _parameters.parallel.planeComm,  MPI_STATUS_IGNORE);
    }
  }
}
