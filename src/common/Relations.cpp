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

#include "Relations.h"
#include "Experiment.h"
#include <iostream>

using namespace boost;
using namespace std;

Relation::Relation()
{
	type[REL_IF] = false;
	type[REL_SS] = false;
	type[REL_SR] = false;
	type[REL_RS] = false;
	type[REL_RR] = false;
}

void Relation::printSingle(bool val, const char* s, fstream& log)
{
	if(!val) {
		log << "__";
	}
	else {
		log << s;
	}
}

void Relation::print(fstream& log)
{
	printSingle(type[REL_IF],"IF",log);
	log << " ";
	printSingle(type[REL_SS],"SS",log);
	log << " ";
	printSingle(type[REL_SR],"SR",log);
	log << " ";
	printSingle(type[REL_RS],"RS",log);
	log << " ";
	printSingle(type[REL_RR],"RR",log);
}

void RelationSet::insert(enum RelType type, int affected, int source)
{
	Link alink = route->getLinkById(affected);
	Link slink = route->getLinkById(source);

	if(alink.source == slink.source) {
		return;
	}

	std::pair<int,int> rel = make_pair(affected,source);

	std::map<std::pair<int,int>,Relation>::iterator i = relations.find(rel);
	if(i == relations.end()) {
		Relation rt;
		rt.affected = affected;
		rt.source = source;
		relations[rel] = rt;
		i = relations.find(rel);
	} 

	i->second.type[type] = true;
}

void RelationSet::insertSet(enum RelType type, int affected, enum Direction dir, vector<int>& set)
{
	for(vector<int>::iterator i = set.begin(); i != set.end(); i++) {
		int e = *i;
		if(e != 0) { // the gateway is no anchor
			if(dir == UP) {
				insert(type,affected,route->getLinkByAnchor(e, true).id);
			}
			else {
				insert(type,affected,route->getLinkByAnchor(e, false).id);
			}
		}
	}
}

void RelationSet::print()
{
	log << relations.size() << std::endl;
	for(std::map<std::pair<int,int>,Relation>::iterator i = relations.begin(); 
			i != relations.end(); i++) {
		log << i->first.first+1 << " " << i->first.second+1 << " ";
		i->second.print(log);
		log << std::endl;
	}
	log.flush();
}

int RelationSet::getRelationsCount()
{
	return relations.size();
}

std::map<std::pair<int,int>,Relation>::iterator RelationSet::begin()
{
	return relations.begin();
}

std::map<std::pair<int,int>,Relation>::iterator RelationSet::end()
{
	return relations.end();
}
