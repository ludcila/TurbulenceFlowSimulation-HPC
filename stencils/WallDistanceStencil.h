#ifndef _WALL_DISTANCE_STENCIL_H_
#define _WALL_DISTANCE_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"


class WallDistanceStencil : public FieldStencil<TurbulentFlowField> {

	protected:
	
		FLOAT _lengthX;
		FLOAT _lengthY;
		FLOAT _lengthZ;

    public:

        /** Constructor
         *
         * @param prefix String with the prefix of the name of the VTK files
         */
        WallDistanceStencil ( const Parameters & parameters );
		~WallDistanceStencil ();

        /** 2D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         */
        virtual void apply ( TurbulentFlowField & flowField, int i, int j );

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        virtual void apply ( TurbulentFlowField & flowField, int i, int j, int k );


};

#endif
