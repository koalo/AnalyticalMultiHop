/*
 * Class for calculating the analytical model
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

#ifndef CALCULATION_H
#define CALCULATION_H

#include <petsc.h>
#include <petscdmnetwork.h>

#include "Queue.h"

struct _p_NODEDATA{
	PetscInt	id;
	PetscScalar   	packet_generation;
	PetscScalar	tx_slots;
	Queue		queue;
};

typedef struct _p_NODEDATA *NODEDATA;

struct _p_SLOTDATA{
  PetscInt	pos;
  PetscScalar   Pactive;
};

typedef struct _p_SLOTDATA *SLOTDATA;

enum {
  VAR_PACTIVE,
  VAR_NVARS
};

typedef struct {
	PetscInt  a;
	PetscInt      outerCircle;
	PetscScalar   inverse;
	TDMASchedule* schedule;
} UserCtx;

typedef struct {
	PetscScalar Paccept;
} Result;

PetscErrorCode EvaluateNode(DM& circuitdm, int v, DMNetworkComponentGenericDataType *arr, const PetscScalar *xarr, UserCtx *user, Result *r, bool debug, PetscScalar* delay = nullptr, PetscScalar* Qmean = nullptr);
PetscErrorCode FormFunction(SNES snes,Vec X, Vec F,void *appctx);
PetscErrorCode SetInitialValues(DM circuitdm,Vec X,void *appctx);
PetscErrorCode FormatResult(Vec& X, DM& circuitdm, UserCtx *user);

#endif

