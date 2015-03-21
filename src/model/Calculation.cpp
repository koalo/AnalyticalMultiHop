/*
 * Class for determining link qualities
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

#include "Calculation.h"
#include <iostream>
#include <fstream>
#include "SqMx.h"

using namespace std;

inline PetscScalar Q(double t, double Q1G) {
	return 1 - pow(1 - Q1G, t);
}

inline PetscScalar geo(double val, double exp) {
	return (1-pow(val,exp))/(1-val);
}

inline PetscScalar PU(PetscScalar a, PetscScalar b) {
	return a + b - a*b;
}

inline PetscScalar PU(PetscScalar a, PetscScalar b, PetscScalar c) {
	return a + b + c - a*b - a*c - b*c + a*b*c;
}

inline PetscScalar PU(PetscScalar a, PetscScalar b, PetscScalar c, PetscScalar d) {
	return a + b + c + d - a*b - a*c - a*d - b*c - b*d - c*d
		+ a*b*c + a*b*d + a*c*d + b*c*d - a*b*c*d;
}

inline PetscScalar PU(PetscScalar a, PetscScalar b, PetscScalar c, PetscScalar d, PetscScalar e) {
	return a + b + c + d + e - a*b - a*c - a*d - a*e - b*c - b*d - b*e - c*d - c*e - d*e
		+ a*b*c + a*b*d + a*b*e + a*c*d + a*c*e + a*d*e + b*c*d  + b*c*e + b*d*e + c*d*e
		- a*b*c*d - a*b*c*e - a*b*d*e - a*c*d*e - b*c*d*e + a*b*c*d*e;
}

#undef __FUNCT__
#define __FUNCT__ "EvaluateLink"
inline PetscErrorCode EvaluateLink(DM& circuitdm, int v, DMCircuitComponentGenericDataType *arr, const PetscScalar *xarr, UserCtx *user, Result *r, bool debug)
{
	PetscInt m = user->m;
	PetscInt mb = user->mb;
	PetscInt m0 = user->m0;
	PetscInt n = user->n;
	PetscScalar L = user->L;
	PetscScalar Lc = user->Lc;
	PetscScalar Ls = user->Ls;
	PetscScalar Lack = user->Lack;

	PetscInt keyv;
	PetscInt offsetlink;
	PetscErrorCode ierr;
	ierr = DMCircuitGetComponentTypeOffset(circuitdm,v,0,&keyv,&offsetlink);CHKERRQ(ierr);

	LINKDATA link = (LINKDATA)(arr+offsetlink);

	PetscScalar PER = link->PER;
	PetscScalar AER = link->AER;

	PetscScalar inflow = 0;
	PetscScalar PnoNeighbourSending = 1;
	PetscScalar PnoHiddenSending = 1;
	PetscScalar PhearsNoACK = 1;
	PetscScalar PXP[8];
	PetscScalar PXA[2];
	PetscScalar PInvBothCollide = 1;
	PetscScalar PInvBothSameCollide = 1;

	for(int i = 0; i < 8; i++) {
		PXP[i] = 1;
	}

	for(int i = 0; i < 2; i++) {
		PXA[i] = 1;
	}

	PetscInt nconnedges;
	const PetscInt *connedges;
	ierr = DMCircuitGetSupportingEdges(circuitdm,v,&nconnedges,&connedges);CHKERRQ(ierr);
	for (PetscInt i = 0; i < nconnedges; i++) {
		PetscInt e = connedges[i];

		const PetscInt *cone;
		ierr = DMCircuitGetConnectedNodes(circuitdm,e,&cone);CHKERRQ(ierr);
		PetscInt vaffected,vsource;
		vaffected = cone[0];
		vsource   = cone[1];

		if (vaffected == v) {
			PetscInt keye;
			PetscInt offsetrel;
			ierr = DMCircuitGetComponentTypeOffset(circuitdm,e,0,&keye,&offsetrel);CHKERRQ(ierr);
			RELATIONDATA relation = (RELATIONDATA)(arr+offsetrel);

			PetscInt offsetsource,offsetsourcelink;
			ierr = DMCircuitGetVariableOffset(circuitdm,vsource,&offsetsource);CHKERRQ(ierr);

			PetscInt keysourcev;
			ierr = DMCircuitGetComponentTypeOffset(circuitdm,vsource,0,
					&keysourcev,&offsetsourcelink);CHKERRQ(ierr);
			LINKDATA sourcelink = (LINKDATA)(arr+offsetsourcelink);

			PetscScalar PnotSending = 1 + xarr[offsetsource+VAR_TAU]*(xarr[offsetsource+VAR_ALPHA] - 1);
			PetscScalar PnotSendingACK = PnotSending;

			if (relation->type[REL_IF]) {
				inflow += xarr[offsetsource+VAR_R]*xarr[offsetsource+VAR_LAMBDA]
					* (1 - sourcelink->arrival_factor);
			}

			bool SS = relation->type[REL_SS];
			bool RS = relation->type[REL_RS];
			bool SR = relation->type[REL_SR];
			bool RR = relation->type[REL_RR];

			if (SS) {
				PnoNeighbourSending *= PnotSending;
			}

			if (RS && !SS) {
				PnoHiddenSending *= PnotSending;
			}

			if (SR) {
				PhearsNoACK *= PnotSendingACK;
			}

			if (SS && RS) {
				if(SR) {
					PInvBothSameCollide *= PnotSending;
				}
			}

			if (RS && !SS) {
				if(SR) {
					PInvBothCollide *= PnotSending;
				}
			}

			if(RS && !SS) {
				PXP[0] *= PnotSending;
			}

			if(SS && RS) {
				PXP[1] *= PnotSending;
			}

			if(RR && !RS && !SR && !SS) {
				PXP[2] *= PnotSendingACK;
			}

			if(RS && RR && !SS && !SR) {
				PXP[3] *= PnotSendingACK;
			}

			if(SS && RR && !SR) {
				PXP[4] *= PnotSendingACK;
			}

			if(SR && RR && !SS && RS) {
				PXP[5] *= PnotSendingACK;
			}

			if(SR && RR && !SS && !RS) {
				PXP[7] *= PnotSendingACK;
			}

			if(SR && RR && SS) {
				PXP[6] *= PnotSendingACK;
			}

			if(SS && !RS) {
				PXA[0] *= PnotSending;
			}

			if(SS && RS) {
				PXA[1] *= PnotSending;
			}
		}
	}

	PetscScalar W0 = pow(2,m0);

	for(int i = 0; i < 8; i++) {
		PXP[i] = 1 - PXP[i];
	}

	for(int i = 0; i < 2; i++) {
		PXA[i] = 1 - PXA[i];
	}

	PXP[0] = Q(2*L,PXP[0]);
	PXP[1] = Q(2,PXP[1]);
	PXP[2] = Q(L+Lack,PXP[2]);
	PXP[3] = Q(Lack+1,PXP[3]);
	PXP[4] = Q(Lack,PXP[4]);
	PXP[5] = Q(2,PXP[5]);
	PXP[6] = Q(1,PXP[6]);
	PXP[7] = Q(2,PXP[7]);
	PXA[0] = Q(Lack,PXA[0]);
	PXA[1] = Q(1,PXA[1]);

	PetscScalar Pcoll = 0;
	PetscScalar PACKlost = 0;

	PetscScalar PCB_l2 = Q(2*L+2,(1 - PInvBothCollide));
	PetscScalar PCB_l1 = 1 - PInvBothSameCollide;

	PetscScalar alphapkt = Q(L,(1-PnoNeighbourSending));
	PetscScalar alphaack = Q(Lack,(1 - PhearsNoACK));

	if(user->handleACKs) {
		// handle acknowledgements
		for(int i = 0; i < 8; i++) {
			Pcoll = Pcoll + PXP[i] - Pcoll*PXP[i];
		}

		for(int i = 0; i < 2; i++) {
			PACKlost = PACKlost + PXA[i] - PACKlost*PXA[i];
		}
	} else {
		// do not handle acknowledgements
		Pcoll = 0;
		for(int i = 0; i <= 1; i++) {
			Pcoll = Pcoll + PXP[i] - Pcoll*PXP[i];
		}

		PACKlost = 0;
		alphaack = 0;
	}

	PetscScalar alpha = alphapkt+alphaack - alphapkt*alphaack;

	Pcoll = Pcoll + (1-Pcoll)*PER;
	PetscScalar PnoACK = Pcoll + (1-Pcoll)*PACKlost;
	PnoACK = PnoACK + (1-PnoACK)*AER;

	link->pcoll = Pcoll;
	link->pnoack = PnoACK;

	PetscScalar y = (PnoACK)*(1-pow(alpha,m));

	PetscScalar packet_generation = link->packet_generation;
	if(user->inverse) {
		PetscInt vStart,vEnd,offset;
		ierr = DMCircuitGetVertexRange(circuitdm,&vStart,&vEnd);CHKERRQ(ierr);
		ierr = DMCircuitGetVariableOffset(circuitdm,vEnd-1,&offset);CHKERRQ(ierr);

		// for avoiding zero pivot
		if(link->packet_generation > 1e-10) {
			packet_generation = xarr[offset];
		}
		else {
			packet_generation = xarr[offset]*1e-9;
		}
	}

	PetscScalar lambda = packet_generation + link->input_factor*inflow;
	PetscScalar q = 1 - exp(-lambda);
	PetscInt md = mb - m0;

	PetscScalar bikj = W0*(1-pow(2*alpha,min(m,md+1)))/(1-2*alpha) + (1-pow(alpha,min(m,md+1)))/(1-alpha) + (pow(2,mb)+1)*pow(alpha,md+1)*(1-pow(alpha,max(0,m-md-1)))/(1-alpha);
	PetscScalar b000 = 1/( 0.5*bikj*geo(y,n+1) + (1-pow(alpha,m))*geo(y,n+1)*(Ls*(1-PnoACK)+Lc*PnoACK) + (1/q)*((1-PnoACK)*(1-pow(alpha,m))*geo(y,n+1) + pow(y,n+1) + pow(alpha,m)*geo(y,n+1)) );

	SqMx<6,PetscScalar> PcrMat = SqMx<6,PetscScalar>::zeros();

	PetscScalar PCR_l2[2];
	PetscScalar PCR_l1[2];

	PCR_l2[0] = 0;
	PCR_l1[0] = 0;
	PCR_l2[1] = user->PCR_l2_1;
	PCR_l1[1] = user->PCR_l1_1;

	PcrMat.a[0][0] = 1;
	PcrMat.a[1][1] = 1;

	PetscScalar am = pow(alpha,m);

	for(int i = 0; i <= 1; i++) {
		for(int j = 0; j <= 1; j++) {
			int s1 = i+j*2+2;
			PcrMat.a[s1][0] = am;
			PcrMat.a[s1][1] = (1-am)*(1-PU(Pcoll,PCR_l2[i],PCR_l1[j]));
		}
	}

	for(int i1 = 0; i1 <= 1; i1++) {
		for(int j1 = 0; j1 <= 1; j1++) {
			for(int i2 = 0; i2 <= 1; i2++) {
				for(int j2 = 0; j2 <= 1; j2++) {
					int s1 = i1+j1*2+2;
					int s2 = i2+j2*2+2;

					if(i2 == 0 && j2 == 0) {
						PcrMat.a[s1][s2] = (1-am)*(Pcoll-PU(PCB_l2,PCB_l1))*(1-PU(PCR_l2[i1],PCR_l1[j1]));
					}
					else if (i2 == 1 || j2 == 1) {
						PetscScalar again_l2 = PU(PCB_l2,PCR_l2[i1]);
						PetscScalar again_l1 = PU(PCB_l1,PCR_l1[j1]);
						PetscScalar result = (1-am);

						if(i2 == 0) {
							result *= 1 - again_l2;
						}
						else {
							result *= again_l2;
						}

						if(j2 == 0) {
							result *= 1 - again_l1;
						}
						else {
							result *= again_l1;
						}

						PcrMat.a[s1][s2] = result;
					}

				}
			}
		}
	}

	SqMx<6,double> PcrMatRetr = PcrMat^(n+1); 

	for(int i = 0; i < 8; i++) {
		r->PXP[i] = PXP[i];
	}

	for(int i = 0; i < 2; i++) {
		r->PXA[i] = PXA[i];
	}

	r->lambda = lambda;

	if(user->betterRetrans) {
		r->Rel = PcrMatRetr.a[2][1];
	}
	else {
		PetscScalar sum = 0;
		PetscScalar am = pow(alpha,m);
		for(int j = 0; j <= n; j++) {
			sum += pow(Pcoll*(1-am),j);
		} 
		r->Rel = 1 - pow(Pcoll*(1-am),n+1) - am*sum;
	}

	r->tau = geo(alpha,m)*geo(y,n+1)*b000;
	r->alphapkt = alphapkt;
	r->alphaack = alphaack;
	r->alpha = alpha;
	r->Pcoll = Pcoll;
	r->PnoACK = PnoACK;
	r->packet_generation = packet_generation;
	link->curR = r->Rel;

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
	PetscInt      offset;
	DMCircuitComponentGenericDataType *arr;

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

	ierr = DMCircuitGetVertexRange(circuitdm,&vStart,&vEnd);CHKERRQ(ierr);
	if(user->inverse) {
		vEnd--;
	}

	ierr = DMCircuitGetComponentDataArray(circuitdm,&arr);CHKERRQ(ierr);

	for (v=vStart; v < vEnd; v++) {
		PetscBool   ghostvtex;
		ierr = DMCircuitIsGhostVertex(circuitdm,v,&ghostvtex);CHKERRQ(ierr);
		if (ghostvtex) {
			continue;
		}

		Result result;
		EvaluateLink(circuitdm,v,arr,xarr,user,&result,false);

		ierr = DMCircuitGetVariableOffset(circuitdm,v,&offset);CHKERRQ(ierr);
		farr[offset+VAR_LAMBDA]     = xarr[offset+VAR_LAMBDA] - result.lambda;
		farr[offset+VAR_R]    = xarr[offset+VAR_R] - result.Rel;
		farr[offset+VAR_TAU]  = xarr[offset+VAR_TAU] - result.tau;
		farr[offset+VAR_ALPHA]  = xarr[offset+VAR_ALPHA] - result.alpha;
	}

	if(user->inverse) {
		// route
		PetscInt vz = (vEnd+vStart)/2-1;

		PetscScalar R = 0;
		PetscInt no = 0;

		for(int z = vz+1; z > vz+1-user->outerCircle; z--) {
			PetscInt v = z;
			PetscBool   ghostvtex;
			ierr = DMCircuitIsGhostVertex(circuitdm,v,&ghostvtex);CHKERRQ(ierr);
			PetscInt to;
			PetscScalar Rx = 1;

			do {
				PetscInt keyv;
				PetscInt offsetlink;
				ierr = DMCircuitGetComponentTypeOffset(circuitdm,v,0,&keyv,&offsetlink);CHKERRQ(ierr);

				LINKDATA link = (LINKDATA)(arr+offsetlink);

				ierr = DMCircuitGetVariableOffset(circuitdm,v,&offset);CHKERRQ(ierr);
				Rx *= xarr[offset+VAR_R]-farr[offset+VAR_R]; // gives result.Rel, see above

				to = link->to;

				v = to+vStart-1;
			} while(to != 0);

			no++;
			R += Rx;
		}

		R /= no;

		ierr = DMCircuitGetVariableOffset(circuitdm,vEnd,&offset);CHKERRQ(ierr);

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
	PetscInt       v, vStart, vEnd, offset;
	Vec            localX;
	PetscScalar    *xarr;
	DMCircuitComponentGenericDataType *arr;
	UserCtx       *user=(UserCtx*)appctx;

	PetscFunctionBegin;
	ierr = DMCircuitGetVertexRange(circuitdm,&vStart,&vEnd);CHKERRQ(ierr);
	if(user->inverse) {
		vEnd--;
	}

	ierr = DMGetLocalVector(circuitdm,&localX);CHKERRQ(ierr);

	ierr = VecSet(X,0.0);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(circuitdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(circuitdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

	PetscScalar firstPGen = 0;
	ierr = DMCircuitGetComponentDataArray(circuitdm,&arr);CHKERRQ(ierr);
	ierr = VecGetArray(localX,&xarr);CHKERRQ(ierr);
	for (v = vStart; v < vEnd; v++) {
		PetscBool   ghostvtex;
		ierr = DMCircuitIsGhostVertex(circuitdm,v,&ghostvtex);CHKERRQ(ierr);
		if (ghostvtex) {
			continue;
		}

		PetscInt offsetlink, offset;
		PetscInt keyv;
		ierr = DMCircuitGetVariableOffset(circuitdm,v,&offset);CHKERRQ(ierr);
		ierr = DMCircuitGetComponentTypeOffset(circuitdm,v,0,&keyv,&offsetlink);CHKERRQ(ierr);

		LINKDATA link = (LINKDATA)(arr+offsetlink);

		xarr[offset+VAR_LAMBDA] = link->packet_generation;
		xarr[offset+VAR_R] = 1;
		xarr[offset+VAR_TAU] = link->packet_generation;
		xarr[offset+VAR_ALPHA] = link->packet_generation;

		if(v == vStart) {
			firstPGen = link->packet_generation;
		}
	}

	if(user->inverse) {
		// extra variable for inverse
		ierr = DMCircuitGetVariableOffset(circuitdm,vEnd,&offset);CHKERRQ(ierr);
		xarr[offset] = firstPGen;
	}

	ierr = VecRestoreArray(localX,&xarr);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(circuitdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(circuitdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(circuitdm,&localX);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

