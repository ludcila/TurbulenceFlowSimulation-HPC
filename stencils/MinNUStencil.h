#ifndef _MIN_NU_STENCIL_H_
#define _MIN_NU_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../TurbulentFlowField.h"


/** this class computes the maximum value of max(velocity)/meshsize for all grid cells.
 *  Originally, one would compute the max. velocity only and adapt it with the meshsize afterwards.
 *  This, however, becomes inconsistent when dealing with non-uniform, e.g. stretched, meshes, since
 *  the meshsize may be different for every grid cell. We therefore determine the max(velocity)/meshsize
 *  and synchronise this value over whole computational domain.
 *  @author Philipp Neumann
 */
class MinNUStencil : public FieldStencil<TurbulentFlowField>, public BoundaryStencil<TurbulentFlowField> {

    private:

        FLOAT _minValue;  //! Stores the maximum module of every component

        /** Sets the maximum value arrays to the value of the cell if it surpasses the current one.
         *
         * 2D version of the function
         * @param flowField Flow field
         * @param i Position in the X direction.
         * @param j Position in the Y direction.
         */
        void cellMinValue(TurbulentFlowField & flowField, int i, int j);

        /** Sets the maximum value arrays to the value of the cell if it surpasses the current one.
         *
         * 3D version of the function
         * @param flowField Flow field
         * @param i Position in the X direction.
         * @param j Position in the Y direction.
         * @param k Position in the Z direction.
         */
        void cellMinValue(TurbulentFlowField & flowField, int i, int j, int k);

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        MinNUStencil (const Parameters & parameters);

        //@ brief Body iterations
        //@{
        void apply (TurbulentFlowField & flowField, int i, int j);
        void apply (TurbulentFlowField & flowField, int i, int j, int k);
        //@}

        //@ brief Boundary iterations for the 2D problem
        //@param flowField Flow field with the state of the fluid
        //@param i Position in the X direction
        //@param j Position in the Y direction
        //@{
        void applyLeftWall   ( TurbulentFlowField & flowField, int i, int j );
        void applyRightWall  ( TurbulentFlowField & flowField, int i, int j );
        void applyBottomWall ( TurbulentFlowField & flowField, int i, int j );
        void applyTopWall    ( TurbulentFlowField & flowField, int i, int j );
        //@}

        //@ brief Boundary iterations for the 3D problem
        //@param flowField Flow field with the state of the fluid
        //@param i Position in the X direction
        //@param j Position in the Y direction
        //@param k Position in the Z direction
        //@{
        void applyLeftWall   ( TurbulentFlowField & flowField, int i, int j, int k );
        void applyRightWall  ( TurbulentFlowField & flowField, int i, int j, int k );
        void applyBottomWall ( TurbulentFlowField & flowField, int i, int j, int k );
        void applyTopWall    ( TurbulentFlowField & flowField, int i, int j, int k );
        void applyFrontWall  ( TurbulentFlowField & flowField, int i, int j, int k );
        void applyBackWall   ( TurbulentFlowField & flowField, int i, int j, int k );
        //@}

        /** Resets the maximum values to zero before computing the timestep
         */
        void reset ();

        /** Returns the array with the maximum modules of the components of the velocity,
         *  divided by the respective local meshsize
         */
        const FLOAT  getMinValue() const;
};

#endif
