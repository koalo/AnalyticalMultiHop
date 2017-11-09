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

#include <petsc.h>
#include <petscdmnetwork.h>

#include "Experiment.h"
#include "Calculation.h"
#include "ResultWriter.h"
#include "Queue.h"

using namespace std;
using namespace boost::filesystem;
using namespace boost::system;

static char help[] = "Run as ./tdma_model --experiment <filename>\n";

int main(int argc, char** argv)
{
	PetscErrorCode       ierr;
	UserCtx user;
	bool debug = true;

	NODEDATA nodes;
	SLOTDATA slots;
	int                  *edges = NULL;

	PetscInt             i;  
	DM                   circuitdm;
	PetscInt             componentkey[2];
	PetscLogStage        stage1,stage2;
	PetscInt             size;

	PetscInt numVertices = 0;
	PetscInt numEdges = 0;

	Experiment experiment;
	char experiment_file[PETSC_MAX_PATH_LEN];

	/* Initialize */
	PetscInitialize(&argc,&argv,NULL,help);

	PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

	/* Create an empty circuit object */
	ierr = DMNetworkCreate(PETSC_COMM_WORLD,&circuitdm);CHKERRQ(ierr);

	/* Register the components in the circuit */
	ierr = DMNetworkRegisterComponent(circuitdm,"linkstruct",
			sizeof(struct _p_NODEDATA),&componentkey[0]);CHKERRQ(ierr);
	ierr = DMNetworkRegisterComponent(circuitdm,"slotsstruct",
			sizeof(struct _p_SLOTDATA),&componentkey[1]);CHKERRQ(ierr);

	ierr = PetscLogStageRegister("Read Data",&stage1);CHKERRQ(ierr);
	PetscLogStagePush(stage1);

	/* READ THE DATA */
	if (!rank) {
		/* Only rank 0 reads the data */

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

		user.outerCircle = experiment.getIntermediate<int>("nodesOnOuterCircle");
		user.inverse = experiment.getParameter<double>("inverse");

		PetscReal freqUp = 1/experiment.getParameter<double>("intervalUp"); // Hz
		freqUp /= 1000000; // uHz

		PetscReal freqDown = 1/experiment.getParameter<double>("intervalDown"); // Hz
		freqDown /= 1000000; // uHz


		/* Read nodes */
		ierr = PetscMalloc1(experiment.getRoute().getNodeCount(), 
				&nodes);CHKERRQ(ierr);

		TDMASchedule& schedule = experiment.getTDMASchedule();
		numVertices = schedule.getNodeCount();
		cout << numVertices << endl;
		for(unsigned int i = 0; i < schedule.getNodeCount(); i++) {
			nodes[i].id = i;

			int slotDuration = experiment.getIntermediate<int>("slotDuration_us");

			// upstream
			if(i > 0) {
				nodes[i].packet_generation = freqUp*slotDuration;
			}
			else {
				nodes[i].packet_generation = 0;
			}

			// TODO no support for downstream for TDMA yet
			
			new (&nodes[i].queue) Queue(); // call constructor
			int K = experiment.getParameter<int>("K");
			nodes[i].queue.create(schedule.getNodes()[i],K);
		}

		/* Read slots */
		int slotsCount = schedule.getTotalTXSlots();
		numEdges = slotsCount;

		// slots for storing mainly the type
		ierr = PetscMalloc1(slotsCount,&slots);CHKERRQ(ierr);
		
		// edges for storing the relation itself
		// *2 for both vertices of an edge
		ierr = PetscMalloc1(2*slotsCount,&edges);CHKERRQ(ierr);

		cout << schedule.getNodes()[0].slots.size() << endl;

		int i = 0;
		for(unsigned int n = 0; n < schedule.getNodeCount(); n++) {
			TDMASchedule::Node& node = schedule.getNodes().at(n);
			nodes[n].tx_slots = 0;
			for(unsigned int pos = 0; pos < node.slots.size(); pos++) {
				auto& slot = node.slots[pos];
				if(slot.type == TDMASchedule::Type::TX) {
					slots[i].pos = pos;
					slots[i].Pactive = 0;

					edges[2*i] = n;
					edges[2*i+1] = slot.counterpart;
					i++;

					nodes[n].tx_slots++;
				}
			}
		}
	}

	user.schedule = &experiment.getTDMASchedule();

	if(debug) {
		cout << "Nodes:" << endl;
		for(int i = 0; i < numVertices; i++) {
			cout << i+1 << endl;
  			cout << nodes[i].packet_generation << endl;
			cout << endl;
		}
		cout << endl;

		cout << "Edges:" << endl;
		for(int i = 0; i < numEdges; i++) {
			cout << edges[2*i] << " " << edges[2*i+1] << endl;
		}
		cout << endl;
	}

	PetscLogStagePop();
	MPI_Barrier(PETSC_COMM_WORLD);

	ierr = PetscLogStageRegister("Create circuit",&stage2);CHKERRQ(ierr);
	PetscLogStagePush(stage2);

	/* Set number of nodes/edges */
	ierr = DMNetworkSetSizes(circuitdm,numVertices,numEdges,
			PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);

	/* Add edge connectivity */
	ierr = DMNetworkSetEdgeList(circuitdm,edges);CHKERRQ(ierr);

	/* Set up the circuit layout */
	ierr = DMNetworkLayoutSetUp(circuitdm);CHKERRQ(ierr);

	if (!rank) {
		ierr = PetscFree(edges);CHKERRQ(ierr);
	}

	/* Add circuit components */
	PetscInt eStart, eEnd, vStart, vEnd;

	ierr = DMNetworkGetVertexRange(circuitdm,&vStart,&vEnd);CHKERRQ(ierr);

	for (i = vStart; i < vEnd; i++) {
		ierr = DMNetworkAddComponent(circuitdm,i,componentkey[0],
				&nodes[i-vStart]);CHKERRQ(ierr);
	}

	ierr = DMNetworkGetEdgeRange(circuitdm,&eStart,&eEnd);CHKERRQ(ierr);
	for (i = eStart; i < eEnd; i++) {
		ierr = DMNetworkAddComponent(circuitdm,i,componentkey[1],
				&slots[i-eStart]);CHKERRQ(ierr);

		/* Add number of variables */
		ierr = DMNetworkAddNumVariables(circuitdm,i,VAR_NVARS);CHKERRQ(ierr);
	}

	if(user.inverse) {
		/* Add additional variable */
		ierr = DMNetworkAddNumVariables(circuitdm,eEnd,1);CHKERRQ(ierr);
	}

	/* Set up DM for use */
	ierr = DMSetUp(circuitdm);CHKERRQ(ierr);

	if (!rank) {
		/* PETSc memcpys the data within DMSetUp. */
	        /* Therefore we can free the data here. */
		ierr = PetscFree(nodes);CHKERRQ(ierr);
		ierr = PetscFree(slots);CHKERRQ(ierr);
	}

	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size > 1) {
		/* Circuit partitioning and distribution of data */
		ierr = DMNetworkDistribute(&circuitdm,0);CHKERRQ(ierr);
	}

	PetscLogStagePop();

	Vec X,F;
	ierr = DMCreateGlobalVector(circuitdm,&X);CHKERRQ(ierr);
	ierr = VecDuplicate(X,&F);CHKERRQ(ierr);

	ierr = SetInitialValues(circuitdm,X,&user);CHKERRQ(ierr);

	SNES snes;

	/* HOOK UP SOLVER */
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
	ierr = SNESSetDM(snes,circuitdm);CHKERRQ(ierr);
	ierr = SNESSetFunction(snes,F,FormFunction,&user);CHKERRQ(ierr);
	ierr = SNESSetType(snes, SNESNRICHARDSON); // faster for our problem, especially by preventing the approximation of the Jacobian
	ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

	ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);

	/* Print out result */
	ResultWriter resultWriter;
	resultWriter.store(X,circuitdm,&user,experiment);

	if(!rank) {
		resultWriter.write(experiment.getResultFileName(experiment_file));
	}

	ierr = VecDestroy(&X);CHKERRQ(ierr);
	ierr = VecDestroy(&F);CHKERRQ(ierr);

	ierr = SNESDestroy(&snes);CHKERRQ(ierr);
	ierr = DMDestroy(&circuitdm);CHKERRQ(ierr);

	PetscFinalize();
	return 0;
}

