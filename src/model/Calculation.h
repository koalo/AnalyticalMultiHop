/*
 * Class for calculating the analytical model
 *
 * Author:	Florian Meier <florian.meier@koalo.de>
 *		Copyright 2015
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
#include <petscdmcircuit.h>

#include "Relations.h"

struct _p_LINKDATA{
  PetscInt	from;
  PetscInt	to;
  PetscScalar   packet_generation;
  PetscScalar   input_factor;
  PetscScalar   arrival_factor;
  PetscScalar pcoll;
  PetscScalar pnoack;
  PetscScalar curR;
  PetscScalar  PER;
  PetscScalar  AER;
};

typedef struct _p_LINKDATA *LINKDATA;

struct _p_RELATIONDATA{
  PetscInt      source;
  PetscInt      affected;
  PetscInt  type[REL_NTYPES];
};

typedef struct _p_RELATIONDATA *RELATIONDATA;

enum {
  VAR_LAMBDA,
  VAR_R,
  VAR_TAU,
  VAR_ALPHA, 
  VAR_NVARS
};

typedef struct {
	PetscInt  m;
	PetscInt  m0;
	PetscInt  mb;
	PetscInt  n;
	PetscScalar  L;
	PetscScalar  Lc;
	PetscScalar  Ls;
	PetscScalar  Lack;
	PetscScalar  TSCinSbs;
	PetscInt     handleACKs;
	PetscInt     betterRetrans;
	PetscScalar     inverse;
	PetscScalar  PCR_l2_1;
	PetscScalar  PCR_l1_1;
	PetscScalar  Sb;
	PetscInt     outerCircle;
} UserCtx;

typedef struct {
	PetscScalar lambda;
	PetscScalar Rel;
	PetscScalar tau;
	PetscScalar alphapkt;
	PetscScalar alphaack;
	PetscScalar alpha;
	PetscScalar Pcoll;
	PetscScalar PXP[8];
	PetscScalar PXA[2];
	PetscScalar PnoACK;
	PetscScalar packet_generation;
} Result;

PetscErrorCode EvaluateLink(DM& circuitdm, int v, DMCircuitComponentGenericDataType *arr, const PetscScalar *xarr, UserCtx *user, Result *r, bool debug);
PetscErrorCode FormFunction(SNES snes,Vec X, Vec F,void *appctx);
PetscErrorCode SetInitialValues(DM circuitdm,Vec X,void *appctx);
PetscErrorCode FormatResult(Vec& X, DM& circuitdm, UserCtx *user);

#endif
