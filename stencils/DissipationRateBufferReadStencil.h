#ifndef __DISSIPATIONRATE_BUFFER_READ_STENCIL_H__
#define __DISSIPATIONEATE_BUFFER_READ_STENCIL_H__

#include "../Stencil.h"
#include "../TurbulentFlowField.h"

class DissipationRateBufferReadStencil : public BoundaryStencil<TurbulentFlowField> {

	private:
	
	FLOAT *_bufferLeftWall;
	FLOAT *_bufferRightWall;
	FLOAT *_bufferTopWall;
	FLOAT *_bufferBottomWall;
	FLOAT *_bufferFrontWall;
	FLOAT *_bufferBackWall;
	
	public:

	DissipationRateBufferReadStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall);
	
	void applyLeftWall(TurbulentFlowField & TurbulentFlowField, int i, int j);

	void applyRightWall(TurbulentFlowField & TurbulentFlowField, int i, int j);

	void applyBottomWall(TurbulentFlowField & TurbulentFlowField, int i, int j);

	void applyTopWall(TurbulentFlowField & TurbulentFlowField, int i, int j);

	void applyLeftWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k);

	void applyRightWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k);

	void applyBottomWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k);

	void applyTopWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k);

	void applyFrontWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k);

	void applyBackWall(TurbulentFlowField & TurbulentFlowField, int i, int j, int k);
	
};


#endif
