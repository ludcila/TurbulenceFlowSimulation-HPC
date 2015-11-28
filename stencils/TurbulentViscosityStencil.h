#ifndef _TURBVISC_STENCIL_H_
#define _TURBVISC_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../TurbulentFlowField.h"



class TurbulentViscosityStencil : public FieldStencil<TurbulentFlowField> {

	

    public:

        /** Constructor
         *
         * @param prefix String with the prefix of the name of the VTK files
         */
        TurbulentViscosityStencil ( const Parameters & parameters );
		~TurbulentViscosityStencil ();

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
