/*
 * Class for calculating a TDMA queue
 *
 * Author:	Florian Kauer <florian.kauer@koalo.de>
 *		Copyright 2017
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Queue.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <petscksp.h>

using namespace std;

void Queue::create(TDMASchedule::Node& nodeSchedule, uint16_t K) {
	this->K = K;
	nodeSchedule.printSlots();
	schedule_length = 0;
	nT = 0;
	for(auto& slot : nodeSchedule.slots) {
		if(slot.type == TDMASchedule::Type::TX) {
			nT++;
			assert(schedule_length < MAX_SCHEDULE_LENGTH);
			T[schedule_length] = schedule_length;
		}
		schedule_length++;
	}

	assert(K <= MAX_K);

	int i = 0;
	uint16_t currentTxId = nT-1;
	for(auto& slot : nodeSchedule.slots) {
		prevTxId[i] = currentTxId; // if the current slot is TX, prevTxId is still the ID of the previous TX slot
		if(slot.type == TDMASchedule::Type::TX) {
			currentTxId = (currentTxId + 1) % nT;
		}
		i++;
	}

	numStates = 0;
	for(uint16_t q = 0; q <= K; q++) {
		for(uint16_t i = 0; i < schedule_length; i++) {
			assert(numStates < MAX_STATES);
			bool tx = nodeSchedule.slots[i].type == TDMASchedule::Type::TX;
			states[numStates] = State{.q=q, .i=i, .tx=tx};
			numStates++;
		}
	}
}

PetscScalar poisson(PetscInt k, PetscScalar l) {
	if(k > 20) {
		return 0;
	}

	uint64_t fac = 1;
	for(int i = 1; i <= k; ++i) {
		fac *= i;
	}
	return (pow(l,k)/fac)*exp(-l);
}

PetscScalar poisson_at_least(PetscInt k, PetscScalar l) {
	PetscScalar p = 1;
	for(int i = 0; i < k; i++) {
		p -= poisson(i,l);
	}
	return p;
}

PetscScalar Ppush(PetscInt k, PetscScalar l, PetscScalar Precv) {
	if(k == 0) {
		return poisson(0,l)*(1-Precv);
	}
	else {
		return Precv*poisson(k-1,l)+(1-Precv)*poisson(k,l);
	}
}

PetscScalar Ppush_at_least(PetscInt k, PetscScalar l, PetscScalar Precv) {
	if(k == 0) {
		return 1;
	}
	else {
		return Precv*poisson_at_least(k-1,l)+(1-Precv)*poisson_at_least(k,l);
	}
}

inline PetscScalar rateaccept(PetscScalar q, PetscInt K, PetscScalar l, PetscScalar Precv) {
	PetscScalar result = Ppush_at_least(K-q,l,Precv)*(K-q);
	for(PetscInt k = 0; k <= K-q-1; k++) {
		result += Ppush(k,l,Precv) * k;
	}
	return result;
}

PetscScalar Queue::calcDelay(PetscInt q, PetscInt i) {
	PetscScalar fq = ceil((q/(PetscScalar)nT)-1);

	PetscScalar result = fq*schedule_length;
	uint16_t packets_left = q; // subtracting fq*nT is not neccessary due to the % nT
	uint16_t nxtTx = prevTxId[i]; 
	nxtTx = (nxtTx+packets_left)%nT;
	uint16_t endi = T[nxtTx];
	if(endi >= i) {
		result += endi - i;
	}
	else {
		result += schedule_length-i+endi;
	}

	result += 1; // transmission duration itself

	return result;
}

PetscErrorCode Queue::addTransition(int idx, uint16_t q, uint16_t i, PetscScalar p) {
	PetscInt j = getPosition(q, i);

	// transpose on the fly, swap i and j
	return MatSetValues(*A,1,&j,1,&idx,&p,ADD_VALUES);
}

PetscErrorCode Queue::calculate(PetscScalar* Paccept, PetscScalar* delay, vector<PetscScalar>* queuedistribution) {
	PetscErrorCode ierr;

	// Create vectors
	Vec x;
	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE,numStates);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);CHKERRQ(ierr);

	Vec b;
	ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

	// Create matrix
	Mat A;
	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,numStates,numStates);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	ierr = MatSetUp(A);CHKERRQ(ierr);

	// Assemble matrix
	setMatrix(&A);

	for(PetscInt i = 0; i < numStates; i++) {
		Queue::State& state = getStates()[i];

		uint8_t v = 0;
		if(state.tx) {
			v = 1;
		}

		PetscScalar Precv = getPrecv(state.i);

		for(int k = 0; k < K-state.q; k++) {
			ierr = addTransition(i,max(state.q-v,0)+k, (state.i+1)%schedule_length,
					-Ppush(k,getLambda(),Precv)); CHKERRQ(ierr);
		}
		ierr = addTransition(i,max(state.q-v,0)+K-state.q,  (state.i+1)%schedule_length,
				-Ppush_at_least(K-state.q,getLambda(),Precv)); CHKERRQ(ierr);

		PetscScalar one = 1.0;
		MatSetValues(A,1,&i,1,&i,&one,ADD_VALUES);
	}

	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = VecSet(b,0);CHKERRQ(ierr);


	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Create the linear solver and set various options
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	KSP ksp;
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp,tolerance,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

	// Set runtime options
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

	// Non-zero guess (required!)
	PetscScalar p = .5;
	ierr = VecSet(x,p);CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Solve the linear system
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	// Setup GMRES with SOR
	PC pc;
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCSOR);CHKERRQ(ierr);
	ierr = KSPSetType(ksp,KSPGMRES);CHKERRQ(ierr);

	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

	KSPConvergedReason reason;
	ierr = KSPGetConvergedReason(ksp, &reason);CHKERRQ(ierr);
	if(reason <= 0) {
		// returns from function
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_CONV_FAILED,"KSP did not converge");
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Handle solution
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	// Normalization
	PetscScalar sum;
	ierr = VecSum(x,&sum);CHKERRQ(ierr);
	ierr = VecScale(x,1/sum);CHKERRQ(ierr);

	PetscScalar* steady;
	ierr = VecGetArray(x,&steady);CHKERRQ(ierr);

	// Calculate Paccept
	*Paccept = 0;
	for(PetscInt i = 0; i < numStates; i++) {
		Queue::State& state = getStates()[i];
		PetscScalar c = steady[i];
		PetscScalar Precv = getPrecv(state.i);
		*Paccept += c * rateaccept(state.q,K,lambda,Precv);
	}
	*Paccept /= getQmean();

	// Calculate delay
	if(delay != nullptr) {
		*delay = 0;
		for(PetscInt i = 0; i < numStates; i++) {
			Queue::State& state = getStates()[i];
			PetscScalar c = steady[i];

			uint8_t v = 0;
			if(state.tx) {
				v = 1;
			}
			(*delay) += c*calcDelay(max(state.q-v,0)+1,(state.i+1)%schedule_length);
		}
	}

	// Calculate queue distribution
	if(queuedistribution != nullptr) {
		queuedistribution->clear();
		queuedistribution->resize(K+1,0);
		for(PetscInt i = 0; i < numStates; i++) {
			Queue::State& state = getStates()[i];
			(*queuedistribution)[state.q] += steady[i];
		}
	}

	// Calculate transmission probabilities
	for(int i = 0; i < schedule_length; i++) {
		// the idle states are at the beginning - so getStates()[i]
		// is the idle state for time step i
		Queue::State& state = getStates()[i];
		if(state.tx) {
			Ptx[i] = 1-steady[i]*schedule_length;
		}
		else {
			Ptx[i] = 0;
		}
	}

	// Free work space
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&b);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

	return 0;
}
