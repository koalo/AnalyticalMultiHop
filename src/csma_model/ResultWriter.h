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

#ifndef RESULTWRITER_H
#define RESULTWRITER_H

#include <petsc.h>
#include <string>
#include <boost/property_tree/ptree.hpp>

#include "Calculation.h"
#include "Route.h"

class ResultWriter {
public:
	ResultWriter();

	PetscErrorCode store(Vec& X, Vec& F, DM& circuitdm, UserCtx *user, Route& route);

	void write(const std::string& filename);

private:
	boost::property_tree::ptree results;
};

#endif

