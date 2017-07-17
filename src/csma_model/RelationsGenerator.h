/*
 * Class for creating the relations for the model
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

#ifndef RELATIONSGENERATOR_H
#define RELATIONSGENERATOR_H

class Experiment;
class Route;
class RelationSet;

#include <vector>

class RelationsGenerator {
	public:
		static void create(Experiment& experiment, Route& route, RelationSet& relations);

	private:
		static void makeUnique(std::vector<int>& v);
};

#endif

