/*
 * Class for writing the results of the analytical model
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

#include "ResultWriter.h"
#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/property_tree/json_parser.hpp>

using namespace std;
using namespace boost;

ResultWriter::ResultWriter() {
}

PetscErrorCode ResultWriter::store(Vec& X, DM& circuitdm, UserCtx *user)
{
	const PetscScalar *xarr;
	PetscErrorCode ierr;
	ierr = VecGetArrayRead(X,&xarr);CHKERRQ(ierr);

	PetscInt vStart, vEnd;
	ierr = DMCircuitGetVertexRange(circuitdm,&vStart,&vEnd);CHKERRQ(ierr);
	if(user->inverse) {
		vEnd--;
	}

	DMCircuitComponentGenericDataType *arr;
	ierr = DMCircuitGetComponentDataArray(circuitdm,&arr);CHKERRQ(ierr);

	int j = 0;
	for (PetscInt v=vStart; v < vEnd; v++, j++) {
		Result result;
		EvaluateLink(circuitdm,v,arr,xarr,user,&result,true);

		property_tree::ptree rtree;
		stringstream strstr;

		rtree.put("lambda",result.lambda);
		rtree.put("Rel",result.Rel);
		rtree.put("tau",result.tau);
		rtree.put("alphapkt",result.alphapkt);
		rtree.put("alphaack",result.alphaack);
		rtree.put("alpha",result.alpha);
		rtree.put("Pcoll",result.Pcoll);
		rtree.put("PnoACK",result.PnoACK);

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

		strstr.str("");
		strstr << j+1;
		results.add_child(strstr.str(),rtree);
	}

	ierr = VecRestoreArrayRead(X,&xarr);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

void ResultWriter::write(const std::string& filename) {
	property_tree::json_parser::write_json(filename,results);
}

