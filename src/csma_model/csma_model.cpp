/*
 * Program for calculating the equation system of
 * the analytical model.
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

#include <iostream>
#include <set>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <petsc.h>
#include <petscdmnetwork.h>

#include "Experiment.h"
#include "Calculation.h"
#include "RelationsGenerator.h"
#include "ResultWriter.h"

using namespace std;
using namespace boost::filesystem;
using namespace boost::system;

static char help[] = "Run as ./csma_model --experiment <filename>\n";

int main(int argc, char** argv)
{
	PetscErrorCode       ierr;
	UserCtx user;
	bool debug = false;

	LINKDATA links;
	RELATIONDATA relations;
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
			sizeof(struct _p_LINKDATA),&componentkey[0]);CHKERRQ(ierr);
	ierr = DMNetworkRegisterComponent(circuitdm,"relationstruct",
			sizeof(struct _p_RELATIONDATA),&componentkey[1]);CHKERRQ(ierr);

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
		PetscInt nodes = experiment.getParameter<int>("nodes");

		user.outerCircle = experiment.getIntermediate<int>("nodesOnOuterCircle");

		int bytesPerACK = 11;
		user.m = experiment.getParameter<int>("MaxNumberOfBackoffs");
		user.m0 = experiment.getParameter<int>("InitialBackoffExponent");
		user.mb = experiment.getParameter<int>("MaxBackoffExponent");
		user.n = experiment.getParameter<int>("MaxNumberOfRetransmissions");

		int broadcast = experiment.getParameter<int>("broadcast");
		if(broadcast && user.n != 0) {
			std::cout << "For Broadcasts, retransmissions are not possible" << std::endl;
			std::cout << help << std::endl;
			return 1;
		}

		PetscScalar LinSym = 8.0*experiment.getParameter<int>("L")/4.0;
		user.L = LinSym/20;
		PetscScalar LackinSym = (bytesPerACK*8.0/4.0);
		user.Lack = (bytesPerACK/10.0);

		// in Symbols
		PetscScalar tACK = 12;
		PetscScalar IFS = 40;
		PetscScalar tmACK = 54; // macAckWaitDuration

		user.Lc = (LinSym+tmACK)/20.0;
		user.Ls = (LinSym+tACK+LackinSym+IFS)/20.0;

		user.handleACKs = experiment.getParameter<int>("handleACKs");

		if(broadcast && user.handleACKs) {
			std::cout << "For Broadcasts, ACKs are not possible" << std::endl;
			std::cout << help << std::endl;
			return 1;
		}

		user.betterRetrans = experiment.getParameter<int>("betterRetrans");
		user.inverse = experiment.getParameter<double>("inverse");
		PetscScalar aUnitBackoffPeriod = 20;
		PetscScalar symDur = 16;
		user.Sb = aUnitBackoffPeriod*symDur;
		PetscScalar phyCCADuration = 8/symDur;
		user.TSCinSbs = (symDur*phyCCADuration)/user.Sb;
		PetscScalar W0 = pow(2,user.m0);

		PetscScalar diff = W0 - user.L;
		if(diff < 0) {
			user.PCR_l2_1 = 1;
		}
		else {
			user.PCR_l2_1 = 1 - ((diff + pow(diff,2))/(W0*W0));
		}
		user.PCR_l1_1 = 1 / W0;

		PetscReal freqUp = 1/experiment.getParameter<double>("intervalUp"); // Hz
		freqUp /= 1000000; // uHz

		PetscReal freqDown = 1/experiment.getParameter<double>("intervalDown"); // Hz
		freqDown /= 1000000; // uHz


		/* Read links */
		ierr = PetscMalloc1(experiment.getRoute().getLinkCount(), 
				&links);CHKERRQ(ierr);

		if(nodes != experiment.getRoute().getLinkCount()/2+1) {
			std::cout << "Number of nodes does not match number of links/2+1" 
				  << std::endl;
			return -1;
		}

		Route& route = experiment.getRoute();
		numVertices = route.getLinkCount();
		for(int i = 0; i < route.getLinkCount(); i++) {
			Link l = route.getLinkById(i);

			links[i].from = l.source;
			links[i].to = l.destination;
  			links[i].pcoll = 0;
  			links[i].pnoack = 0;
  			links[i].curR = 0;

			double BER = route.getBER(i);
			links[i].PER = 1 - pow(1 - BER,experiment.getParameter<int>("L")*8);
			links[i].AER = (1-pow(1 - BER,bytesPerACK*8));

			if(!broadcast) {
				if(l.up) {
					// upstream
					links[i].input_factor = 1;
					links[i].arrival_factor = 0;
				}
				else {
					// downstream
					int descendantsChild = route.getDescendants(l.destination);
					int descendantsParent = route.getDescendants(l.source);
					links[i].input_factor = (1+descendantsChild) /
								((double)descendantsParent);
					links[i].arrival_factor = 1.0/(1.0+descendantsChild);
				}

				if(l.up) {
					// upstream
					links[i].packet_generation = freqUp*user.Sb;
				}
				else if(links[i].from == 0) {
					// downstream - only inner edges
					links[i].packet_generation = freqDown * user.Sb *
								links[i].input_factor*(nodes-1);
				}
				else {
					links[i].packet_generation = 0;
				}
			}
			else {
				// For broadcasts, the data is equally distributed over all downstreams
				// and does not spread!
				links[i].input_factor = 0;
				links[i].arrival_factor = 1;

				if(l.up) {
					links[i].packet_generation = 0;
				}
				else {
					int children = route.getDirectChildren(l.source);
					links[i].packet_generation = freqDown*user.Sb/children;
					cout << l.source << " " << children << " " << links[i].packet_generation << endl;
				}
			}
		}

		/* Generate Relations */
		RelationsGenerator::create(experiment, experiment.getRoute(),
				experiment.getRelations());

		/* Read relations */
		RelationSet relationSet = experiment.getRelations();
		int relationsCount = relationSet.getRelationsCount();
		numEdges = relationsCount;
		if(user.inverse) {
			// edge from every link to the node handling the inverse variable
			numEdges += route.getLinkCount();
		}

		// relations for storing mainly the type
		ierr = PetscMalloc1(relationsCount,&relations);CHKERRQ(ierr);
		
		// edges for storing the relation itself
		// *2 for both vertices of an edge
		ierr = PetscMalloc1(2*numEdges,&edges);CHKERRQ(ierr);

		int i = 0;
		for(auto& rel : relationSet) {
			relations[i].affected = rel.second.affected;
			relations[i].source = rel.second.source;

			for(int t = 0; t < REL_NTYPES; t++) {
				relations[i].type[t] = rel.second.type[t];
			}

			edges[2*i] = relations[i].affected;
			edges[2*i+1] = relations[i].source;
			i++;
		}

		assert(i == relationsCount);

		if(user.inverse) {
			PetscInt links = route.getLinkCount();
			for(int i = 0; i < links; i++) {
				edges[2*relationsCount+2*i] = links; // id of inverse variable node
				edges[2*relationsCount+2*i+1] = i;
			}
		}
	}

	if(debug) {
		cout << "Links:" << endl;
		for(int i = 0; i < numVertices; i++) {
			cout << i << endl;
  			cout << links[i].from << endl;
  			cout << links[i].to << endl;
  			cout << links[i].packet_generation << endl;
  			cout << links[i].input_factor << endl;
  			cout << links[i].arrival_factor << endl;
  			cout << links[i].PER << endl;
  			cout << links[i].AER << endl;
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
	if(user.inverse) {
		numVertices++;  // for extra variable (for inverse) 
	}

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
	if(user.inverse) {
		vEnd--;
	}

	for (i = vStart; i < vEnd; i++) {
		ierr = DMNetworkAddComponent(circuitdm,i,componentkey[0],
				&links[i-vStart]);CHKERRQ(ierr);

		/* Add number of variables */
		ierr = DMNetworkAddNumVariables(circuitdm,i,VAR_NVARS);CHKERRQ(ierr);
	}

	if(user.inverse) {
		// extra variable for inverse
		ierr = DMNetworkAddNumVariables(circuitdm,vEnd,1);CHKERRQ(ierr);
	}

	ierr = DMNetworkGetEdgeRange(circuitdm,&eStart,&eEnd);CHKERRQ(ierr);
	if(user.inverse) {
		Route& route = experiment.getRoute();
		eEnd -= route.getLinkCount();
	}

	for (i = eStart; i < eEnd; i++) {
		ierr = DMNetworkAddComponent(circuitdm,i,componentkey[1],
				&relations[i-eStart]);CHKERRQ(ierr);
	}

	/* Set up DM for use */
	ierr = DMSetUp(circuitdm);CHKERRQ(ierr);

	if (!rank) {
		/* PETSc memcpys the data within DMSetUp. */
	        /* Therefore we can free the data here. */
		ierr = PetscFree(links);CHKERRQ(ierr);
		ierr = PetscFree(relations);CHKERRQ(ierr);
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
	ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

	ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);

	/* Print out result */
	ResultWriter resultWriter;
	Route& route = experiment.getRoute();
	resultWriter.store(X,F,circuitdm,&user,route);

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

