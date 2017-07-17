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

#ifndef RELATIONS_H
#define RELATIONS_H

#include <map>
#include <vector>

class Experiment;
class Route;
class RelationsGenerator;

enum RelType {
	REL_IF,
	REL_SS,
	REL_SR,
	REL_RS,
	REL_RR,
	REL_NTYPES
};

class Relation {
	public:
		int affected;
		int source;
		bool type[REL_NTYPES];

		Relation();
		void printSingle(bool val, const char* s);
		void print();
};

class RelationSet {
	public:
		int getRelationsCount();

		std::map<std::pair<int,int>,Relation>::iterator begin();
		std::map<std::pair<int,int>,Relation>::iterator end();

		void print();

	private:
		enum Direction {
			UP,
			DOWN
		};

		Route* route;

		int nodes;
		std::map<char,int> typeCount;
		std::map<std::pair<int,int>,Relation> relations;

		void insert(enum RelType type, int affected, int source, bool checkMatch = true);
		void insertSet(enum RelType type, int affected, enum Direction dir, std::vector<int>& set);

		friend RelationsGenerator;
};

#endif

