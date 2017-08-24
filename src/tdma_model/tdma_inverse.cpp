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

static char help[] = "Run as ./tdma_inverse --queue QUEUE_LENGTH --tx 0,1,0,0,...,1 (--Paccept TARGET_PACCEPT | --delay TARGET_DELAY)\n\n";

#include <petscsnes.h>
#include <iostream>
#include <iomanip>
#include "Queue.h"

using namespace std;

extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode FormInitialGuess(Vec);

struct Context {
	Queue queue;
	PetscScalar PacceptTarget;
	PetscScalar delayTarget;
	PetscScalar delayTargetNormalized;
	PetscScalar delayMin;
	PetscScalar delayMax;
};

int main(int argc,char **argv)
{
	SNES           snes;                   /* SNES context */
	Mat		 J;
	Vec            x,r;             /* vectors */
	PetscErrorCode ierr;
	PetscInt       its;
	PetscMPIInt    size;

	ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"This is a uniprocessor example only!");

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Read parameters 
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	Context context;
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
	string s = buffer;

	TDMASchedule::Node schedule;
	schedule.TXfromCommaSeparatedString(s);

	PetscBool foundPaccept;
	ierr = PetscOptionsGetScalar(PETSC_NULL,NULL,"--Paccept",
			&context.PacceptTarget,&foundPaccept);CHKERRQ(ierr);

	PetscBool foundDelay;
	ierr = PetscOptionsGetScalar(PETSC_NULL,NULL,"--delay",
			&context.delayTarget,&foundDelay);CHKERRQ(ierr);

	if((foundDelay && foundPaccept) || !(foundDelay || foundPaccept)) {
		std::cout << "Either provide Paccept OR delay!" << std::endl;
		std::cout << help << std::endl;
		return 1;
	}

	if(foundPaccept) {
		context.delayTarget = -1;
	}
	else {
		context.PacceptTarget = -1;
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Setup schedule 
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	context.queue.setTolerance(1e-10);
	context.queue.create(schedule,queueLength);
	for(int i = 0; i < context.queue.getScheduleLength(); i++) {
		context.queue.setPrecv(i,0);
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Create nonlinear solver context
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Calculate delay boundaries
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if(context.delayTarget > 0) {
		PetscScalar Paccept;
		context.queue.setLambda(0);
		ierr = context.queue.calculate(&Paccept,&context.delayMin); CHKERRQ(ierr);
		context.queue.setLambda(1e10);
		ierr = context.queue.calculate(&Paccept,&context.delayMax); CHKERRQ(ierr);

		if(context.delayTarget > context.delayMax) {
			cout << "Target delay larger than maximum delay " << context.delayMax << "!" << endl;
			cout << "Lambda: inf Paccept: nan Delay: nan" << endl;
			return 0;
		}
		else if(context.delayTarget < context.delayMin) {
			cout << "Target delay smaller than minimum delay " << context.delayMin << "!" << endl;
			return 1;
		}

		context.delayTargetNormalized = (context.delayTarget-context.delayMin)/(context.delayMax-context.delayMin);

		ierr = SNESSetTolerances(snes,0.001,1e-8,1e-8,50,10000); CHKERRQ(ierr);
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Create vector data structures; set function evaluation routine
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	   Note that we form 1 vector from scratch and then duplicate as needed.
	 */
	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE,1);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&r);CHKERRQ(ierr);


	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Create matrix data structure; set Jacobian evaluation routine
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	ierr = MatCreate(PETSC_COMM_WORLD,&J);CHKERRQ(ierr);
	ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,1,1);CHKERRQ(ierr);
	ierr = MatSetFromOptions(J);CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(J,3,NULL);CHKERRQ(ierr);

	/*
	   Set function evaluation routine and vector
	 */
	ierr = SNESSetFunction(snes,r,FormFunction,&context);CHKERRQ(ierr);
	ierr = SNESSetJacobian(snes,J,J,SNESComputeJacobianDefault,(void*)FormFunction);CHKERRQ(ierr);
	ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);


	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Evaluate initial guess; then solve nonlinear system
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = FormInitialGuess(x);CHKERRQ(ierr);
	ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);
	ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SNES iterations = %D\n\n",its);CHKERRQ(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Get result
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	SNESConvergedReason reason;
	ierr = SNESGetConvergedReason(snes, &reason);CHKERRQ(ierr);
	if(reason > 0) {
		const PetscScalar *xx;
		ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
		PetscScalar lambda = exp(xx[0]);
		PetscScalar Paccept;
		PetscScalar delay;

		context.queue.setLambda(lambda);
		ierr = context.queue.calculate(&Paccept,&delay); CHKERRQ(ierr);

		cout << "Lambda: " << exp(xx[0]) << " Paccept: " << Paccept << " Delay: " << delay << endl;
	}
	else {
		cout << "No solution found! " << reason << endl;
	}
	

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Clean up
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&r);CHKERRQ(ierr);
	ierr = SNESDestroy(&snes);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return ierr;
}

PetscErrorCode FormInitialGuess(Vec x)
{
	PetscErrorCode ierr;
	PetscScalar    pfive = .50;
	ierr = VecSet(x,pfive);CHKERRQ(ierr);
	return 0;
}

PetscScalar sigInv(PetscScalar y) {
	return log(y/(1-y));
}

PetscErrorCode FormFunction(SNES snes,Vec x,Vec f,void *ctx)
{
	Context* context = (Context*)ctx;
	const PetscScalar *xx;
	PetscScalar       *ff;
	PetscErrorCode    ierr;

	ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
	ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

	PetscScalar lambda = exp(xx[0]);
	PetscScalar delay;
	PetscScalar Paccept;

	cout << lambda << endl;
	context->queue.setLambda(lambda);
	ierr = context->queue.calculate(&Paccept,&delay); CHKERRQ(ierr);

	if(context->PacceptTarget > 0) {
		double shift = 1e-10; // avoid sigInv(1)
		ff[0] = sigInv(Paccept-shift) - sigInv(context->PacceptTarget-shift);
	}
	else {
		PetscScalar delayNormalized = (delay-context->delayMin)/(context->delayMax-context->delayMin);
		ff[0] = sigInv(delayNormalized) - sigInv(context->delayTargetNormalized);
	}

	cout << setprecision(10) << "xx[0] " << xx[0] << " lambda " << lambda << " Paccept " << Paccept << " delay " << delay << " err " << ff[0] << endl;

	ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
	ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);
	return 0;
}
