#ifndef _TURBVISC_KEPS_STENCIL_H_
#define _TURBVISC_KEPS_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../TurbulentFlowField.h"
#include "MixingLengthModel.h"



class TurbulentViscosityKepsStencil : public FieldStencil<TurbulentFlowField> {

	
    private:

      
          FLOAT _localVelocity           [ 27 * 3 ];
        FLOAT _localMeshsize           [ 27 * 3 ];
        FLOAT _localTurbulentViscosity [ 27 * 3 ]; 
	FLOAT _localK [27*3];
	FLOAT _localeps[27*3];

        MixingLengthModel * _mixingLength;
        
        
    public:

        /** Constructor
         *
         * @param prefix String with the prefix of the name of the VTK files
         */
        TurbulentViscosityKepsStencil ( const Parameters & parameters );
		~TurbulentViscosityKepsStencil ();

        /** 2D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         */
        void apply ( TurbulentFlowField & flowField, int i, int j );

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        void apply ( TurbulentFlowField & flowField, int i, int j, int k );

   
};

#endif
