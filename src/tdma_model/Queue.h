/*
 * Class for calculating a TDMA queue
 *
 * Author:	Florian Kauer <florian.kauer@koalo.de>
 *		Copyright 2015-2017
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

#ifndef QUEUE_H
#define QUEUE_H

#include <vector>
#include <cstdint>
#include <iostream>
#include <map>
#include <cassert>

#include <petscsys.h>
#include <petscmat.h>

#include "TDMASchedule.h"

//#define MAX_SCHEDULE_LENGTH 64
#define MAX_SCHEDULE_LENGTH 128
#define MAX_K 30
#define MAX_STATES (MAX_SCHEDULE_LENGTH*(MAX_K+1))

class Queue {
public:
	class State {
	public:
        uint16_t q,i;
	bool tx;
	};

	uint32_t getPosition(uint16_t q, uint16_t i) {
		return q*schedule_length+i;
	}

	PetscErrorCode calculate(PetscScalar* Paccept, PetscScalar* delay = nullptr, std::vector<PetscScalar>* queuedistribution = nullptr);

	uint16_t getK() const {
		return K;
	}

	PetscErrorCode addTransition(int idx, uint16_t q, uint16_t i, PetscScalar p);

	void setMatrix(Mat* matrix) {
		this->A = matrix;
	}

	void create(TDMASchedule::Node& nodeSchedule, uint16_t K);

	PetscScalar getLambda() {
		return lambda;
	}

	PetscScalar getPrecv(uint16_t i) {
		return Precv[i];
	}

	void setPrecv(uint16_t i, PetscScalar p) {
		Precv[i] = p;
	}

	PetscScalar getPtx(uint16_t i) {
		return Ptx[i];
	}

	void setLambda(PetscScalar l) {
		lambda = l;
	}

	uint16_t getScheduleLength() {
		return schedule_length;
	}

	void setTolerance(double tol) {
		tolerance = tol;
	}

	PetscScalar calcDelay(PetscInt q, PetscInt i);

	PetscScalar getQmean() {
		PetscScalar Qmean = 0;
		for(int i = 0; i < schedule_length; i++) {
			PetscScalar Precv = getPrecv(i);
			Qmean += (1-Precv)*lambda+Precv*(lambda+1);
		}
		Qmean /= schedule_length;
		return Qmean;
	}

private:
	std::array<State,MAX_STATES> states;
	PetscScalar lambda;
	std::array<uint16_t,MAX_SCHEDULE_LENGTH> prevTxId;
	std::array<uint16_t,MAX_SCHEDULE_LENGTH> T;
	std::array<PetscScalar,MAX_SCHEDULE_LENGTH> Ptx;
	std::array<PetscScalar,MAX_SCHEDULE_LENGTH> Precv;
	uint16_t K = 0;
	double tolerance = 1e-3;
	uint16_t schedule_length;
	uint16_t nT;
	uint16_t numStates;
	Mat* A;

	std::array<State,MAX_STATES>& getStates() {
		return states;
	}
};

#endif

