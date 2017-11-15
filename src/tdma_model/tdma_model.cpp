/*
 * Program for calculating the equation system of
 * the analytical TDMA model.
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

#include <iostream>
#include <set>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <petsc.h>

#include "Experiment.h"
#include "Queue.h"
#include "Route.h"

using namespace std;
using namespace boost::filesystem;
using namespace boost::system;

static char help[] = "Run as ./tdma_model --experiment <filename>\n";

class Calculation {
public:
	Calculation(Experiment& experiment)
        : experiment(experiment) {
		inverse = experiment.getParameter<double>("inverse");
	}

	PetscErrorCode calculate(PetscScalar freqUp) {
		PetscErrorCode ierr;
		int slotDuration = experiment.getIntermediate<int>("slotDuration_us");

		/* Build data structures */
		TDMASchedule& schedule = experiment.getTDMASchedule();
		visited.clear();
		visited.resize(schedule.getNodeCount(), false);
		calculated.clear();
		calculated.resize(schedule.getNodeCount(), false);
		pendingVisits.clear();
		pendingVisits.push_back(0);
		pendingCalculations.clear();
		Pactive.clear();
		Pactive.resize(schedule.getNodeCount(),vector<PetscScalar>(schedule.getNodes()[0].slots.size()));
		results.clear();
		results.resize(schedule.getNodeCount());

		/* Run */
		// forward visits
		while(!pendingVisits.empty()) {
			int id = pendingVisits.back();
			pendingVisits.pop_back();
			forward(id,experiment.getRoute());
		}

		// calculate and backtrack
		while(!pendingCalculations.empty()) {
			int id = pendingCalculations.back();
			assert(visited[id] == true);
			assert(calculated[id] == false);
			pendingCalculations.pop_back();

			// upstream
			PetscScalar packet_generation = freqUp*slotDuration;
			// no support for downstream for TDMA yet

			int K = experiment.getParameter<int>("K");
			ierr = calculateNode(id, K, packet_generation, schedule.getNodes()[id], &(results[id])); CHKERRQ(ierr);
			results[id].delay *= slotDuration / 1000000.0; // slotDuration in us
			
			calculated[id] = true;
			
			if(id == 0) {
				assert(pendingCalculations.empty());
			}
			else {
				int p = experiment.getRoute().getPredecessor(id);
				backtrack(p,experiment.getRoute());
			}
		}

		// Evaluate results
		resultPtree.clear();
		for(unsigned int j = 0; j < results.size(); j++) {
			Result& result = results[j];
			boost::property_tree::ptree rtree;
			stringstream strstr;

			rtree.put("Paccept",result.Paccept);
			rtree.put("Qmean",result.qmean);

			// Calculate path reliability and delay
			if(j != 0) {
				int nextHop = j;
				double Rtotal = 1;
				double Dtotal = 0;
				int hops = 0;
				int firstHop = -1;
				do {
					Rtotal *= results[nextHop].Paccept;
					Dtotal += results[nextHop].delay;
					hops++;

					TDMASchedule::Node& node = experiment.getTDMASchedule().getNodes().at(nextHop);
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

				results[j].Rtotal = Rtotal;

				rtree.put("Rtotal",Rtotal);
				rtree.put("D",results[j].delay);
				rtree.put("Dtotal",Dtotal);
				rtree.put("hops",hops);
				rtree.put("perHopDelay",Dtotal/hops);
				rtree.put("firstHop",firstHop);

				PetscScalar freqGen = result.packet_generation/slotDuration; // uHz
				freqGen *= 1000000; // Hz
				rtree.put("intervalGen",1/freqGen);
			}

			// Fin
			strstr.str("");
			strstr << j;
			resultPtree.add_child(strstr.str(),rtree);
		}

		PetscFunctionReturn(0);
	}

	PetscScalar getInverse() {
		return inverse;
	}

	PetscScalar reliabilityForOuterCircle() {
		int outerCircle = experiment.getIntermediate<int>("nodesOnOuterCircle");

		PetscScalar sum = 0;
		for(unsigned int j = results.size()-outerCircle; j < results.size(); j++) {
			sum += results[j].Rtotal;
		}
		return sum/outerCircle;
	}

	void write_results(const std::string& filename) {

		boost::property_tree::json_parser::write_json(filename,resultPtree);
	}

private:
	typedef struct {
		PetscScalar Paccept;
		PetscScalar packet_generation;
		PetscScalar delay;
		PetscScalar qmean;
		PetscScalar Rtotal;
	} Result;

	void forward(int id, Route& route) {
		assert(visited[id] == false);
		assert(calculated[id] == false);
		visited[id] = true;
		
		// push all children
		bool childrenFound = false;
		for(auto n = route.getAdjacentVertices(id); n.first != n.second; n.first++) {
			const auto& c = *(n.first);
			bool incoming = route.getPredecessor(c) == id;
			if(incoming) {
				// is a child
				childrenFound = true;
				assert(visited[c] == false);
				assert(calculated[c] == false);
				pendingVisits.push_back(c);
			}
		}

		if(!childrenFound) {
			// calculation can be executed
			pendingCalculations.push_back(id);
		}
	}

	void backtrack(int id, Route& route) {
		assert(visited[id] == true);
		assert(calculated[id] == false);

		// check if all children were calculated
		for(auto n = route.getAdjacentVertices(id); n.first != n.second; n.first++) {
			const auto& c = *(n.first);
			bool incoming = route.getPredecessor(c) == id;
			if(incoming) {
				// is a child
				if(!calculated[c]) {
					return;
				}
			}
		}

		// all finished, calculation can be executed
		pendingCalculations.push_back(id);
	}

	PetscErrorCode calculateNode(int id, int K, PetscScalar packet_generation, TDMASchedule::Node& node_schedule, Result *r)
	{
		PetscErrorCode ierr;
		Queue queue;
		queue.create(node_schedule,K);

		for(unsigned int pos = 0; pos < node_schedule.slots.size(); pos++) {
			auto& slot = node_schedule.slots[pos];
			if(slot.type == TDMASchedule::Type::RX) {
				queue.setPrecv(pos, Pactive[slot.counterpart][pos]);
			}
			else {
				queue.setPrecv(pos, 0);
			}
		}

		if(id > 0) {
			queue.setLambda(packet_generation);
			r->packet_generation = packet_generation;
			ierr = queue.calculate(&r->Paccept,&r->delay); CHKERRQ(ierr);
		}
		else {
			queue.setLambda(0);
			r->packet_generation = 0;
			r->Paccept = 1; // the sink does not need to queue
			r->delay = 0;
		}

		r->qmean = queue.getQmean();

		for(unsigned int pos = 0; pos < node_schedule.slots.size(); pos++) {
			auto& slot = node_schedule.slots[pos];
			if(slot.type == TDMASchedule::Type::TX) {
				Pactive[id][pos] = queue.getPtx(pos);
			}
		}

		PetscFunctionReturn(0);
	}

	Experiment& experiment;
	vector<bool> visited;
	vector<bool> calculated;
	vector<int> pendingVisits;
	vector<int> pendingCalculations;
	vector<vector<PetscScalar>> Pactive;
	vector<Result> results;
	PetscScalar inverse;
	boost::property_tree::ptree resultPtree;
};

PetscScalar sigInv(PetscScalar y) {
	return log(y/(1-y));
}

PetscErrorCode FormInverseFunction(SNES snes,Vec x,Vec f,void *ctx)
{
	Calculation* calculator = (Calculation*)ctx;
	const PetscScalar *xx;
	PetscScalar       *ff;
	PetscErrorCode    ierr;

	ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
	ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

	PetscScalar freqUp = exp(xx[0]);

	ierr = calculator->calculate(freqUp); CHKERRQ(ierr);

	ff[0] = sigInv(calculator->getInverse()) - sigInv(calculator->reliabilityForOuterCircle());

	ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
	ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

int main(int argc, char** argv)
{
	PetscErrorCode ierr;

	Experiment experiment;
	char experiment_file[PETSC_MAX_PATH_LEN];

	/* Initialize */
	PetscInitialize(&argc,&argv,NULL,help);

	/* Read data */
	PetscBool found;
	ierr = PetscOptionsGetString(PETSC_NULL,NULL,"--experiment",
			experiment_file,PETSC_MAX_PATH_LEN-1,&found);CHKERRQ(ierr);

	if(!found) {
		std::cout << "No experiment file given" << std::endl;
		std::cout << help << std::endl;
		return 1;
	}

	experiment.read(experiment_file);

	if(!experiment.getTDMASchedule().isCalculated()) {
		std::cout << "No TDMA Schedule generated" << std::endl;
		std::cout << help << std::endl;
		return 1;
	}

	PetscScalar freqUp = 1/experiment.getParameter<double>("intervalUp"); // Hz
	freqUp /= 1000000; // uHz

	Calculation calculator(experiment);
	if(!calculator.getInverse()) {
		ierr = calculator.calculate(freqUp); CHKERRQ(ierr);
	}
	else {
		SNES snes;
		Vec x,r;
		Mat J;

		// Create
		ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
		ierr = VecSetSizes(x,PETSC_DECIDE,1);CHKERRQ(ierr);
		ierr = VecSetFromOptions(x);CHKERRQ(ierr);
		ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
		ierr = MatCreate(PETSC_COMM_WORLD,&J);CHKERRQ(ierr);
		ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,1,1);CHKERRQ(ierr);

		// Setup
		ierr = VecSet(x,1);CHKERRQ(ierr); // initial guess
		ierr = MatSetFromOptions(J);CHKERRQ(ierr);
		ierr = MatSeqAIJSetPreallocation(J,3,NULL);CHKERRQ(ierr);
		ierr = SNESSetFunction(snes,r,FormInverseFunction,&calculator);CHKERRQ(ierr);
		ierr = SNESSetJacobian(snes,J,J,SNESComputeJacobianDefault,(void*)FormInverseFunction);CHKERRQ(ierr);
		ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

		// Calculate
		ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);

		// Release
		ierr = VecDestroy(&x);CHKERRQ(ierr);
		ierr = VecDestroy(&r);CHKERRQ(ierr);
		ierr = SNESDestroy(&snes);CHKERRQ(ierr);
	}
	calculator.write_results(experiment.getResultFileName(experiment_file));

	PetscFinalize();
	return 0;
}

