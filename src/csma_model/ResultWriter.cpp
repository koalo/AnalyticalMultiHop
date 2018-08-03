/*
 * Class for writing the results of the analytical model
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

#include "ResultWriter.h"
#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/property_tree/json_parser.hpp>

using namespace std;
using namespace boost;

ResultWriter::ResultWriter() {
}

PetscErrorCode ResultWriter::store(Vec& X, Vec& F, DM& circuitdm, UserCtx *user, Route& route)
{
	const PetscScalar *xarr;
	const PetscScalar *farr;
	PetscErrorCode ierr;
	ierr = VecGetArrayRead(X,&xarr);CHKERRQ(ierr);
	ierr = VecGetArrayRead(F,&farr);CHKERRQ(ierr);

	PetscInt vStart, vEnd;
	ierr = DMNetworkGetVertexRange(circuitdm,&vStart,&vEnd);CHKERRQ(ierr);
	if(user->inverse) {
		vEnd--;
	}

	DMNetworkComponentGenericDataType *arr;
	ierr = DMNetworkGetComponentDataArray(circuitdm,&arr);CHKERRQ(ierr);

	int x = 0;
	vector<PetscScalar> delayVec;
	delayVec.resize(vEnd-vStart);
	vector<property_tree::ptree> rtrees;
	rtrees.resize(vEnd-vStart);
	for (PetscInt v=vStart; v < vEnd; v++, x++) {
		Result result;
		EvaluateLink(circuitdm,v,arr,xarr,user,&result,true);

		property_tree::ptree& rtree = rtrees[v-vStart];
		stringstream strstr;

		Link l = route.getLinkById(x);
		rtree.put("from",l.source);
		rtree.put("to",l.destination);
		rtree.put("up",l.up);

		rtree.put("lambda",result.lambda);
		rtree.put("Rel",result.Rel);
		rtree.put("tau",result.tau);
		rtree.put("alphapkt",result.alphapkt);
		rtree.put("alphaack",result.alphaack);
		rtree.put("alpha",result.alpha);
		rtree.put("Pcoll",result.Pcoll);
		rtree.put("PnoACK",result.PnoACK);
		rtree.put("queueAccept",result.queueAccept);
#ifdef DELAY
		rtree.put("EDsl",result.EDsl);
		rtree.put("EDml",result.EDml);
		rtree.put("EDnl",result.EDnl);
#endif
		rtree.put("cfdrop",result.cfdrop);

		PetscScalar intervalGen = 1/((result.packet_generation/user->Sb)*1000000);
		if(intervalGen > 999999999) {
			intervalGen = 999999999;
		}
		rtree.put("intervalGen", intervalGen);

		for(unsigned int i = 0; i < sizeof(result.PXP)/sizeof(result.PXP[0]); i++) {
			strstr.str("");
			strstr << "PXP[" << i << "]";
			rtree.put(strstr.str(), result.PXP[i]);
		}

		for(unsigned int i = 0; i < sizeof(result.PXA)/sizeof(result.PXA[0]); i++) {
			strstr.str("");
			strstr << "PXA[" << i << "]";
			rtree.put(strstr.str(), result.PXA[i]);
		}


		// Calculate path reliability
		PetscScalar Rtotal;
		ierr = CalculatePathReliability(circuitdm,arr,xarr,farr,Rtotal,v,vStart,(vEnd-vStart)/2,l.up);CHKERRQ(ierr);
		rtree.put("Rtotal",Rtotal);

		// Calculate successful sending delay (does not consider betterRetrans!)
		PetscScalar EDsl = 0;
		for(PetscInt j = 0; j <= user->n; j++) {
			PetscScalar PrCjC = (1-result.Pcoll*(1-pow(result.alpha,user->m+1))) * pow(result.Pcoll,j) * pow((1-pow(result.alpha,user->m)),j);
			PrCjC = PrCjC / (1 - pow(result.Pcoll*(1-pow(result.alpha,user->m+1)),user->n+1));
			PetscScalar EDslj = user->Ls + j*user->Lc;
			for(PetscInt h = 0; h <= j; h++) {
				PetscScalar Eth = user->TSCinSbs;
				for(PetscInt i = 0; i <= user->m; i++) {
					PetscScalar inner = i*user->TSCinSbs;
					for(PetscInt k = 0; k <= i; k++) {
						PetscScalar Wk = 0;
						if(k == 0) {
							Wk = pow(2,user->m0);
						}
						else if(k <= (user->mb - user->m0)) {
							Wk = pow(2,k+user->m0);
						}
						else {
							Wk = pow(2,user->mb);
						}
						inner += ((Wk-1)/2); 
					}
					Eth += (pow(result.alpha,i)*(1-result.alpha)/(1-pow(result.alpha,user->m+1))) * inner;
				}
				EDslj += Eth;
			}
			EDsl += PrCjC * EDslj;
		}

		delayVec[x] = EDsl*user->Sb/1000000.0; // Sb in s
	}

	x = 0;
	for (PetscInt v=vStart; v < vEnd; v++, x++) {
		property_tree::ptree& rtree = rtrees[v-vStart];

		// Calculate path delay
		PetscScalar Dtotal = 0;
		PetscInt hops = 0;
		PetscInt firstHop = -1;

		PetscInt to;
		PetscInt vn = v;

		do {
			PetscInt keyv;
			PetscInt offsetlink;
			ierr = DMNetworkGetComponentKeyOffset(circuitdm,vn,0,&keyv,&offsetlink);CHKERRQ(ierr);

			LINKDATA link = (LINKDATA)(arr+offsetlink);

			Dtotal += delayVec[vn-vStart];
			hops++;

			to = link->to;
			if(firstHop == -1) {
				firstHop = to;
			}

			vn = to+vStart-1;
		} while(to != 0);

		rtree.put("D",delayVec[v-vStart]);
		rtree.put("Dtotal",Dtotal);
		rtree.put("hops",hops);
		rtree.put("perHopDelay",Dtotal/hops);
		rtree.put("firstHop",firstHop);

		stringstream strstr;
		strstr.str("");
		strstr << x+1;
		results.add_child(strstr.str(),rtree);
	}

	ierr = VecRestoreArrayRead(X,&xarr);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

void ResultWriter::write(const std::string& filename) {
	property_tree::json_parser::write_json(filename,results);
}

