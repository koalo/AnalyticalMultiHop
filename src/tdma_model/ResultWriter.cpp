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

PetscErrorCode ResultWriter::store(Vec& X, DM& circuitdm, UserCtx *user, TDMASchedule& schedule, Experiment& experiment)
{
	const PetscScalar *xarr;
	PetscErrorCode ierr;
	ierr = VecGetArrayRead(X,&xarr);CHKERRQ(ierr);

	PetscInt vStart, vEnd;
	ierr = DMNetworkGetVertexRange(circuitdm,&vStart,&vEnd);CHKERRQ(ierr);

	DMNetworkComponentGenericDataType *arr;
	ierr = DMNetworkGetComponentDataArray(circuitdm,&arr);CHKERRQ(ierr);

	vector<Result> resultVec;
	resultVec.resize(vEnd-vStart);
	vector<PetscScalar> delayVec;
	delayVec.resize(vEnd-vStart);
	vector<PetscScalar> qmeanVec;
	qmeanVec.resize(vEnd-vStart);
	int slotDuration = experiment.getIntermediate<int>("slotDuration_us");
	unsigned int j = 0;
	for (PetscInt v=vStart; v < vEnd; v++, j++) {
		EvaluateNode(circuitdm,v,arr,xarr,user,&(resultVec[j]),true,&(delayVec[j]),&(qmeanVec[j]));
		delayVec[j] *= slotDuration / 1000000.0; // slotDuration in us
	}

	for(unsigned int j = 0; j < resultVec.size(); j++) {
		Result& result = resultVec[j];
		property_tree::ptree rtree;
		stringstream strstr;

		rtree.put("Paccept",result.Paccept);
		rtree.put("Qmean",qmeanVec[j]);

		// Calculate path reliability and delay
		if(j != 0) {
			int nextHop = j;
			double Rtotal = 1;
			double Dtotal = 0;
			int hops = 0;
			int firstHop = -1;
			do {
				Rtotal *= resultVec[nextHop].Paccept;
				Dtotal += delayVec[nextHop];
				hops++;

				TDMASchedule::Node& node = schedule.getNodes().at(nextHop);
				nextHop = -1;
				for(unsigned int pos = 0; pos < node.slots.size(); pos++) {
					auto& slot = node.slots[pos];
					if(slot.type == TDMASchedule::Type::TX) {
						assert(nextHop == -1 || nextHop == slot.counterpart);
						nextHop = slot.counterpart;
					}
				}
				if(firstHop == -1) {
					firstHop = nextHop;
				}
				assert(nextHop != -1);
			} while (nextHop != 0);

			rtree.put("Rtotal",Rtotal);
			rtree.put("D",delayVec[j]);
			rtree.put("Dtotal",Dtotal);
			rtree.put("hops",hops);
			rtree.put("perHopDelay",Dtotal/hops);
			rtree.put("firstHop",firstHop);
		}

		// Fin
		strstr.str("");
		strstr << j;
		results.add_child(strstr.str(),rtree);
	}

	ierr = VecRestoreArrayRead(X,&xarr);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

void ResultWriter::write(const std::string& filename) {
	property_tree::json_parser::write_json(filename,results);
}

