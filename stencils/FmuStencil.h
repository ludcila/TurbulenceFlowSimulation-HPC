#ifndef _Fmu_STENCIL_H_
#define _Fmu_STENCIL_H_

#include "../TurbulentFlowField.h"
#include "../Stencil.h"
#include "../Parameters.h"
#include "StencilFunctions.h"

/** Field stencil to compute the right hand side of the pressure equation.
 */
class FmuStencil : public FieldStencil<TurbulentFlowField> {

     private:

        FLOAT _localMeshsize           [ 27 * 3 ];
        FLOAT _localK                  [ 27 * 3 ];
	FLOAT _localeps                [ 27 * 3 ];
    public:
        //constructor
        FmuStencil ( const Parameters & parameters );
        //apply methods
        void apply ( TurbulentFlowField & flowField, int i, int j );
        void apply ( TurbulentFlowField & flowField, int i, int j, int k );
};

#endif
