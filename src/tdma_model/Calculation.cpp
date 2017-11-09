/*
 * Class for calculating nodes
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

#include "Calculation.h"
#include <iostream>
#include <fstream>
#include "Queue.h"

using namespace std;

#undef __FUNCT__
#define __FUNCT__ "EvaluateNode"
PetscErrorCode EvaluateNode(DM& circuitdm, int v, DMNetworkComponentGenericDataType *arr, const PetscScalar *xarr, UserCtx *user, Result *r, bool debug, PetscScalar* delay, PetscScalar* Qmean)
{
	PetscInt keyv;
	PetscInt offsetlink;
	PetscErrorCode ierr;
	ierr = DMNetworkGetComponentKeyOffset(circuitdm,v,0,&keyv,&offsetlink);CHKERRQ(ierr);

	NODEDATA node = (NODEDATA)(arr+offsetlink);

	PetscInt nconnedges;

	if(!user->inverse) {
		node->queue.setLambda(node->packet_generation);
	}
	else {
		PetscInt eStart,eEnd,offset;
		ierr = DMNetworkGetVertexRange(circuitdm,&eStart,&eEnd);CHKERRQ(ierr);
		ierr = DMNetworkGetVariableOffset(circuitdm,eEnd,&offset);CHKERRQ(ierr);

		// for avoiding zero pivot TODO remove
		//if(node->packet_generation > 1e-10) {
		//	node->queue.setLambda(xarr[offset]);
		//}
		//else {
			node->queue.setLambda(xarr[offset]*1e-9);
		//}
	}

	const PetscInt *connedges;
	ierr = DMNetworkGetSupportingEdges(circuitdm,v,&nconnedges,&connedges);CHKERRQ(ierr);
	for (PetscInt i = 0; i < nconnedges; i++) {
		PetscInt e = connedges[i];

		const PetscInt *cone;
		ierr = DMNetworkGetConnectedVertices(circuitdm,e,&cone);CHKERRQ(ierr);
		PetscInt from;
		from = cone[0];

		PetscInt keye;
		PetscInt offsetrel;
		ierr = DMNetworkGetComponentKeyOffset(circuitdm,e,0,&keye,&offsetrel);CHKERRQ(ierr);
		SLOTDATA slot = (SLOTDATA)(arr+offsetrel);

		if(v == from) {
			// TX slot -> will be handled in Queue
		}
		else {
			// RX slot
			PetscInt offsettx,offsettxnode,keytxv;
			ierr = DMNetworkGetVariableOffset(circuitdm,from,&offsettx);CHKERRQ(ierr);
			ierr = DMNetworkGetComponentKeyOffset(circuitdm,from,0,
					&keytxv,&offsettxnode);CHKERRQ(ierr);
			node->queue.setPrecv(slot->pos, slot->Pactive);
		}
	}

	if(node->id > 0) {
		ierr = node->queue.calculate(&r->Paccept,delay); CHKERRQ(ierr);
	}
	else {
		r->Paccept = 1; // the sink does not need to queue
		if(delay) {
			*delay = 0;
		}
	}

	if(Qmean) {
		*Qmean = node->queue.getQmean();
	}

	// write back Pactive
	for (PetscInt i = 0; i < nconnedges; i++) {
		PetscInt e = connedges[i];

		const PetscInt *cone;
		ierr = DMNetworkGetConnectedVertices(circuitdm,e,&cone);CHKERRQ(ierr);
		PetscInt from;
		from = cone[0];

		PetscInt keye;
		PetscInt offsetrel;
		ierr = DMNetworkGetComponentKeyOffset(circuitdm,e,0,&keye,&offsetrel);CHKERRQ(ierr);
		SLOTDATA slot = (SLOTDATA)(arr+offsetrel);

		if(v == from) {
			// TX slot
			slot->Pactive = node->queue.getPtx(slot->pos);
		}
	}

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormFunction"
PetscErrorCode FormFunction(SNES snes,Vec X, Vec F,void *appctx)
{
	PetscErrorCode ierr;
	DM             circuitdm;
	Vec           localX,localF;
	PetscInt      v,vStart,vEnd;

	UserCtx       *user=(UserCtx*)appctx;

	const PetscScalar *xarr;
	PetscScalar   *farr;
	DMNetworkComponentGenericDataType *arr;

	PetscFunctionBegin;
	ierr = SNESGetDM(snes,&circuitdm);CHKERRQ(ierr);
	ierr = DMGetLocalVector(circuitdm,&localX);CHKERRQ(ierr);
	ierr = DMGetLocalVector(circuitdm,&localF);CHKERRQ(ierr);
	ierr = VecSet(F,0.0);CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(circuitdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(circuitdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(circuitdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(circuitdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);

	ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);
	ierr = VecGetArray(localF,&farr);CHKERRQ(ierr);

	ierr = DMNetworkGetVertexRange(circuitdm,&vStart,&vEnd);CHKERRQ(ierr);

	ierr = DMNetworkGetComponentDataArray(circuitdm,&arr);CHKERRQ(ierr);

	vector<Result> resultVec; // only for inverse calculation
	if(user->inverse) {
		resultVec.resize(vEnd-vStart);
	}

	unsigned int j = 0;
	for (v=vStart; v < vEnd; v++, j++) {
		PetscBool   ghostvtex;
		ierr = DMNetworkIsGhostVertex(circuitdm,v,&ghostvtex);CHKERRQ(ierr);
		if (ghostvtex) {
			continue;
		}

		if(!user->inverse) {
			Result result;
			EvaluateNode(circuitdm,v,arr,xarr,user,&result,false);
		}
		else {
			EvaluateNode(circuitdm,v,arr,xarr,user,&(resultVec[j]),false);
		}

		PetscInt nconnedges;
		const PetscInt *connedges;
		ierr = DMNetworkGetSupportingEdges(circuitdm,v,&nconnedges,&connedges);CHKERRQ(ierr);
		for (PetscInt i = 0; i < nconnedges; i++) {
			PetscInt e = connedges[i];

			const PetscInt *cone;
			ierr = DMNetworkGetConnectedVertices(circuitdm,e,&cone);CHKERRQ(ierr);
			PetscInt from;
			from = cone[0];

			PetscInt keye;
			PetscInt offsetrel;
			ierr = DMNetworkGetComponentKeyOffset(circuitdm,e,0,&keye,&offsetrel);CHKERRQ(ierr);
			SLOTDATA slot = (SLOTDATA)(arr+offsetrel);

			if(v == from) {
				PetscInt      offset;
				ierr = DMNetworkGetVariableOffset(circuitdm,e,&offset);CHKERRQ(ierr);
				// TX Slot
				farr[offset+VAR_PACTIVE]     = xarr[offset+VAR_PACTIVE] - slot->Pactive;
			}
		}
	}

	if(user->inverse) {
		PetscScalar R = 0;
		PetscInt no = 0;

		for(unsigned int j = resultVec.size()-user->outerCircle; j < resultVec.size(); j++) {
			PetscScalar Rx;
			int nextHop = j;
			double Rtotal = 1;
			do {
				Rtotal *= resultVec[nextHop].Paccept;

				TDMASchedule::Node& node = user->schedule->getNodes().at(nextHop);
				nextHop = -1;
				for(unsigned int pos = 0; pos < node.slots.size(); pos++) {
					auto& slot = node.slots[pos];
					if(slot.type == TDMASchedule::Type::TX) {
						assert(nextHop == -1 || nextHop == slot.counterpart);
						nextHop = slot.counterpart;
					}
				}
				assert(nextHop != -1);
			} while (nextHop != 0);
			no++;
			R += Rx;
		}

		R /= no;

		PetscInt eStart,eEnd,offset;
		ierr = DMNetworkGetVertexRange(circuitdm,&eStart,&eEnd);CHKERRQ(ierr);
		ierr = DMNetworkGetVariableOffset(circuitdm,eEnd,&offset);CHKERRQ(ierr);

		// in order to avoid zero pivot
		if(R >= 0.9999) {
			R += xarr[offset];
		}
		else {
			R += 0.00001*xarr[offset];
		}

		farr[offset] = R-user->inverse;
	}

	ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
	ierr = VecRestoreArray(localF,&farr);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(circuitdm,&localX);CHKERRQ(ierr);

	ierr = DMLocalToGlobalBegin(circuitdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(circuitdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(circuitdm,&localF);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetInitialValues"
PetscErrorCode SetInitialValues(DM circuitdm,Vec X,void *appctx)
{
	PetscErrorCode ierr;
	PetscInt       eStart, eEnd;
	Vec            localX;
	PetscScalar    *xarr;
	UserCtx       *user=(UserCtx*)appctx;
	DMNetworkComponentGenericDataType *arr;

	PetscFunctionBegin;
	ierr = DMNetworkGetEdgeRange(circuitdm,&eStart,&eEnd);CHKERRQ(ierr);

	ierr = DMGetLocalVector(circuitdm,&localX);CHKERRQ(ierr);

	ierr = VecSet(X,0.0);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(circuitdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(circuitdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

	ierr = DMNetworkGetComponentDataArray(circuitdm,&arr);CHKERRQ(ierr);
	ierr = VecGetArray(localX,&xarr);CHKERRQ(ierr);
	for (PetscInt e = eStart; e < eEnd; e++) {
		PetscInt offset;
		ierr = DMNetworkGetVariableOffset(circuitdm,e,&offset);CHKERRQ(ierr);

		xarr[offset+VAR_PACTIVE] = 0.01;
	}

	if(user->inverse) {
		// extra variable for inverse
		PetscInt eStart,eEnd,offset;
		ierr = DMNetworkGetVertexRange(circuitdm,&eStart,&eEnd);CHKERRQ(ierr);
		ierr = DMNetworkGetVariableOffset(circuitdm,eEnd,&offset);CHKERRQ(ierr);
		xarr[offset] = 10;
	}

	ierr = VecRestoreArray(localX,&xarr);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(circuitdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(circuitdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(circuitdm,&localX);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

