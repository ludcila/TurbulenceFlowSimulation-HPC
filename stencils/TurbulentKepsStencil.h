#ifndef _KEPS_STENCIL_H_
#define _KEPS_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../TurbulentFlowField.h"

/** Stencil to compute the k and eps
 */
class TurbulentKepsStencil : public FieldStencil<TurbulentFlowField> {

	private:
	
		FLOAT _localVelocity           [ 27 * 3 ];
		FLOAT _localMeshsize           [ 27 * 3 ];
		FLOAT _localTurbulentViscosity [ 27 * 3 ]; 
		FLOAT _localK                  [ 27 * 3 ];
		FLOAT _localeps                [ 27 * 3 ];
		FLOAT _localfmu                [ 27 * 3 ];

    public:

        /** Constructor
         * @param parameters Parameters of the problem
         */
        TurbulentKepsStencil(const Parameters & parameters);

        /** Apply the stencil in 2D
         * @param flowField Flow field information
         * @param i Position in the X direction
         * @param j Position in the Y direction
         */
        void apply ( TurbulentFlowField & flowField, int i, int j );

        /** Apply the stencil in 3D
         * @param flowField Flow field information
         * @param i Position in the X direction
         * @param j Position in the Y direction
         * @param k Position in the Z direction
         */
        void apply ( TurbulentFlowField & flowField, int i, int j, int k );
};

#endif
