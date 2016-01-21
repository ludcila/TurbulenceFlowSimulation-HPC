#include "KepsObstacleStencil.h"

KepsObstacleStencil::KepsObstacleStencil(const Parameters & parameters) :
	FieldStencil<TurbulentFlowField> (parameters) {}

void KepsObstacleStencil::apply ( TurbulentFlowField & flowField, int i, int j){

	const int obstacle = flowField.getFlags().getValue(i, j);

	if ((obstacle & OBSTACLE_SELF) == 1){
	
		/* if top cell is fluid */
		if ((obstacle & OBSTACLE_TOP) == 0){
			flowField.getKineticEnergy().getScalar(i, j) = -flowField.getKineticEnergy().getScalar(i, j+1) * _parameters.meshsize->getDy(i, j) / _parameters.meshsize->getDy(i, j+1);
			flowField.getTurbulentViscosity().getScalar(i, j) = -flowField.getTurbulentViscosity().getScalar(i, j+1) * _parameters.meshsize->getDy(i, j) / _parameters.meshsize->getDy(i, j+1);
			flowField.getDissipationRate().getScalar(i, j) = flowField.getDissipationRate().getScalar(i, j+1);
		}
		/* if bottom cell is fluid */
		else if ((obstacle & OBSTACLE_BOTTOM) == 0){
			flowField.getDissipationRate().getScalar(i, j) = flowField.getDissipationRate().getScalar(i, j-1);
		}
		/* if right cell is fluid */
		else if ((obstacle & OBSTACLE_RIGHT) == 0){
			flowField.getKineticEnergy().getScalar(i, j) = -flowField.getKineticEnergy().getScalar(i+1, j)  * _parameters.meshsize->getDx(i, j) / _parameters.meshsize->getDx(i+1, j);
			flowField.getTurbulentViscosity().getScalar(i, j) = -flowField.getTurbulentViscosity().getScalar(i+1, j)  * _parameters.meshsize->getDx(i, j) / _parameters.meshsize->getDx(i+1, j);
			flowField.getDissipationRate().getScalar(i, j) = flowField.getDissipationRate().getScalar(i+1, j);
		}
		/* if left cell is fluid */
		else if ((obstacle & OBSTACLE_LEFT) == 0){
			flowField.getDissipationRate().getScalar(i, j) = flowField.getDissipationRate().getScalar(i-1, j);
		}
		
	}



}

void KepsObstacleStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k){

	const int obstacle = flowField.getFlags().getValue(i, j, k);

	if ((obstacle & OBSTACLE_SELF) == 1){
	
		/* if top cell is fluid */
		if ((obstacle & OBSTACLE_TOP) == 0){
			flowField.getKineticEnergy().getScalar(i, j, k) = 0;
			flowField.getTurbulentViscosity().getScalar(i, j, k) = 0;
			flowField.getDissipationRate().getScalar(i, j, k) = flowField.getDissipationRate().getScalar(i, j+1, k);
		}
		/* if bottom cell is fluid */
		if ((obstacle & OBSTACLE_BOTTOM) == 0){
			flowField.getDissipationRate().getScalar(i, j, k) = flowField.getDissipationRate().getScalar(i, j-1, k);
		}
		/* if right cell is fluid */
		if ((obstacle & OBSTACLE_RIGHT) == 0){
			flowField.getKineticEnergy().getScalar(i, j, k) = 0;
			flowField.getTurbulentViscosity().getScalar(i, j, k) = 0;
			flowField.getDissipationRate().getScalar(i, j, k) = flowField.getDissipationRate().getScalar(i+1, j,k);
		}
		/* if left cell is fluid */
		if ((obstacle & OBSTACLE_LEFT) == 0){
			flowField.getDissipationRate().getScalar(i, j, k) = flowField.getDissipationRate().getScalar(i-1, j, k);
		}
		/* if top and right cells are fluid (top-right corner)*/
		if ((obstacle & OBSTACLE_TOP) == 0 && (obstacle & OBSTACLE_RIGHT) == 0){
			flowField.getDissipationRate().getScalar(i, j, k) = (flowField.getDissipationRate().getScalar(i+1, j, k) + flowField.getDissipationRate().getScalar(i, j+1, k)) * 0.5;
		}
		
	}

}
