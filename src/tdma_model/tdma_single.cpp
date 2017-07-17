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

static char help[] = "Run as ./tdma_single --queue QUEUE_LENGTH --gen 0.5 --tx 0,1,0,....,0 -Precv 0.2,0,...,0.3\n\n";

#include <petscsnes.h>
#include <iostream>
#include <iomanip>
#include "Queue.h"
#include <regex>

using namespace std;

int main(int argc,char **argv)
{
	PetscErrorCode ierr;

	ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Read parameters 
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscInt queueLength;

	PetscBool found;
	ierr = PetscOptionsGetInt(PETSC_NULL,NULL,"--queue",
			&queueLength,&found);CHKERRQ(ierr);
	if(!found) {
		std::cout << "Queue length has to be provided!" << std::endl;
		std::cout << help << std::endl;
		return 1;
	}

	char buffer[500];
	ierr = PetscOptionsGetString(PETSC_NULL,NULL,"--tx",
			buffer,sizeof(buffer),&found);CHKERRQ(ierr);
	if(!found) {
		std::cout << "TX slots have to be provided!" << std::endl;
		std::cout << help << std::endl;
		return 1;
	}


	TDMASchedule::Node schedule;

	string s = buffer;
	regex re("[\\s,]+");
	sregex_token_iterator it(s.begin(), s.end(), re, -1);
	sregex_token_iterator reg_end;
	for (; it != reg_end; ++it) {
		int tx = atoi(it->str().c_str());
		TDMASchedule::Slot slot;
		slot.counterpart = 0;
		if(tx) {
			slot.type = TDMASchedule::Type::TX;
		}
		else {
			slot.type = TDMASchedule::Type::IDLE;
		}
		schedule.slots.push_back(slot);
	}

	Queue queue;
	queue.setTolerance(1e-10);
	queue.create(schedule,queueLength);

	PetscScalar lambda;
	ierr = PetscOptionsGetScalar(PETSC_NULL,NULL,"--gen",&lambda,&found);
	if(!found) {
		std::cout << "gen has to be provided!" << std::endl;
		std::cout << help << std::endl;
		return 1;
	}
	queue.setLambda(lambda);

	ierr = PetscOptionsGetString(PETSC_NULL,NULL,"--Precv",
			buffer,sizeof(buffer),&found);CHKERRQ(ierr);
	if(!found) {
		std::cout << "The Precv have to be provided!" << std::endl;
		std::cout << help << std::endl;
		return 1;
	}

	s = buffer;
	sregex_token_iterator it2(s.begin(), s.end(), re, -1);
	sregex_token_iterator reg_end2;
	unsigned int i = 0;
	for (; it2 != reg_end2; ++it2) {
		double Precv = atof(it2->str().c_str());
		queue.setPrecv(i,Precv);
		i++;
	}
	if(i != schedule.slots.size()) {
		std::cout << "--tx and --Precv have to be of equal length!" << std::endl;
		std::cout << help << std::endl;
		return 1;
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Calculate 
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscScalar Paccept;
	PetscScalar delay;
	vector<PetscScalar> queuedistribution;
	ierr = queue.calculate(&Paccept,&delay,&queuedistribution); CHKERRQ(ierr);

	for(unsigned int q = 0; q < queuedistribution.size(); q++) {
		cout << "Queue " << q << " " << queuedistribution[q] << endl;
	}

	cout << "Paccept " << Paccept << " delay " << delay << endl;

	ierr = PetscFinalize();
	return ierr;
}

