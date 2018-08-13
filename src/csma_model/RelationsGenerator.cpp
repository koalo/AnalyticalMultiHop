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
#include "RelationsGenerator.h"
#include "Experiment.h"

using namespace boost;
using namespace std;

void RelationsGenerator::makeUnique(vector<int>& v) {
	sort(v.begin(),v.end());
	unique(v.begin(),v.end());
}

void RelationsGenerator::create(Experiment& experiment, Route& route, RelationSet& relations)
{
	relations.nodes = experiment.getParameter<int>("nodes");
	relations.route = &route;

	std::vector<int> SSupup;
	std::vector<int> SSdownup;
	std::vector<int> SSupdown;
	std::vector<int> SSdowndown;

	// iterate over all vertices
	// for each vertex handle the link from the direction of the source to this vertex
	// and the link into the direction of the sink 
	auto ai = route.getAllVertices();
	for (; ai.first != ai.second; ai.first++) {
		const graph_traits<Route::Graph>::vertex_descriptor& anchor = *(ai.first);

		SSupup.clear();
		SSdownup.clear();
		SSupdown.clear();
		SSdowndown.clear();

		// the gateway is no anchor
		if(anchor == 0) {
			continue;
		}

		// handle the uplink
		Link uplink = route.getLinkByAnchor(anchor, true);
		auto v2i = route.getAdjacentVertices(uplink.source);
		for(; v2i.first != v2i.second; v2i.first++) {
			const graph_traits<Route::Graph>::vertex_descriptor& v2 = *(v2i.first);

			// v2 is in range of the sender of the uplink (v1)
			// therefore, each uplink with anchor v2 is in SS
			// This does not include two links outgoing of
			// anchor, but there are none in upstream, anyway!
			SSupup.push_back(v2);

			// if two uplinks are cascaded, mark them as inflow edges for each other
			Link uplink2 = route.getLinkByAnchor(v2, true);
			if(uplink.source == uplink2.destination) {
				relations.insert(REL_IF,uplink.id,uplink2.id);
			}

			// search for all w2 in range of v2 so that v2->w2 is a downlink
			// This includes downlinks from anchor, too!
			auto w2i = route.getAdjacentVertices(v2);
			for(; w2i.first != w2i.second; w2i.first++) {
				const graph_traits<Route::Graph>::vertex_descriptor& w2 = *(w2i.first);
				Link downlink2 = route.getLinkByAnchor(w2, false);

				if(downlink2.source == v2) {
					SSupdown.push_back(w2);
				}
			}
		}

		// handle the downlink
		Link downlink = route.getLinkByAnchor(anchor, false);
		v2i = route.getAdjacentVertices(downlink.source);
		for(; v2i.first != v2i.second; v2i.first++) {
			const graph_traits<Route::Graph>::vertex_descriptor& v2 = *(v2i.first);

			// v2 is in range of the sender of the downlink (v1)
			// therefore, each uplink with anchor v2 is in SS
			// This includes anchor, too!
			SSdownup.push_back(v2);

			// search for all w2 in range of v2 so that v2->w2 is a downlink
			auto w2i = route.getAdjacentVertices(v2);
			for(; w2i.first != w2i.second; w2i.first++) {
				const graph_traits<Route::Graph>::vertex_descriptor& w2 = *(w2i.first);
				Link downlink2 = route.getLinkByAnchor(w2, false);

				if(downlink2.source == v2) {
					// Does not include other downlinks
					// from downlink.source!
					SSdowndown.push_back(w2);
				}

				// if two downlinks are cascaded, mark them as inflow edges for each other
				if(downlink2.source == downlink.destination) {
					relations.insert(REL_IF,downlink2.id,downlink.id);
				}
			}

			// downlinks with the same source
			if(v2 != anchor) {
				Link downlink2 = route.getLinkByAnchor(v2, false);
				if(downlink2.source == downlink.source) {
					SSdowndown.push_back(v2);
				}
			}
		}

		makeUnique(SSupup);
		makeUnique(SSdownup);
		makeUnique(SSupdown);
		makeUnique(SSdowndown);

		relations.insertSet(REL_SS,uplink.id,RelationSet::UP,SSupup); 
		relations.insertSet(REL_SS,uplink.id,RelationSet::DOWN,SSupdown); 
		relations.insertSet(REL_SS,downlink.id,RelationSet::UP,SSdownup); 
		relations.insertSet(REL_SS,downlink.id,RelationSet::DOWN,SSdowndown); 

		relations.insertSet(REL_RS,uplink.id,RelationSet::UP,SSdownup); 
		relations.insertSet(REL_RS,uplink.id,RelationSet::DOWN,SSdowndown); 
		relations.insertSet(REL_RS,downlink.id,RelationSet::UP,SSupup); 
		relations.insertSet(REL_RS,downlink.id,RelationSet::DOWN,SSupdown); 

		relations.insertSet(REL_SR,uplink.id,RelationSet::UP,SSupdown); 
		relations.insertSet(REL_SR,uplink.id,RelationSet::DOWN,SSupup); 
		relations.insertSet(REL_SR,downlink.id,RelationSet::UP,SSdowndown); 
		relations.insertSet(REL_SR,downlink.id,RelationSet::DOWN,SSdownup); 

		relations.insertSet(REL_RR,uplink.id,RelationSet::UP,SSdowndown); 
		relations.insertSet(REL_RR,uplink.id,RelationSet::DOWN,SSdownup); 
		relations.insertSet(REL_RR,downlink.id,RelationSet::UP,SSupdown); 
		relations.insertSet(REL_RR,downlink.id,RelationSet::DOWN,SSupup); 
	}

	relations.print();
}

