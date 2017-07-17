/*
 * Class for creating a TDMA schedule
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

#include "TDMAGenerator.h"
#include "Experiment.h"
#include <iostream>

using namespace boost;
using namespace std;

static const int slotDurationDSME = 480*16; // microseconds / slot -> DSME
static const int slotDurationTSCH = 10000; // microseconds / slot -> TSCH (10 ms)

int TDMAGenerator::requiredSlots(Route& route, int node) {
	int desc = route.getDescendants(node);
	if(node == 0) {
		// root only needs RX slots
		return desc;
	}
	else {
		// RX+TX for every desc + 1 TX for own traffic
		return 2*desc+1;
	}
}

class TASCNode {
public:
	virtual void handleNode(int nextSlot) {
		int v = nextChild();
		if(v != -1) {
			markChildVisited();
			nodes->at(v)->forward(nextSlot);
		}
		else if(id != 0) { // node is leaf or sub tree fully handled
			for(int slot = nextSlot; slot <= nextSlot + si - 1; slot++) {
				assign(slot,parent,TDMASchedule::Type::TX,11);
				nodes->at(parent)->assign(slot,id,TDMASchedule::Type::RX,11);
			}
			auto& node = schedule->getNodes()[id];
			node.printSlots();
			nodes->at(parent)->backtrack(nextSlot+si,si);
		}
	}

	void forward(int nextSlot) {
		si = 1;
		handleNode(nextSlot);
	}

	void backtrack(int nextSlot, int slotsOfChild) {
		if(id != 0) {
			si = si + slotsOfChild;
		}
		handleNode(nextSlot);
	}

	virtual void assign(int slot, int counterpart, TDMASchedule::Type type, int channel) {
		schedule->getNodes()[id].slots[slot].type = type;
		schedule->getNodes()[id].slots[slot].counterpart = counterpart;
		schedule->getNodes()[id].slots[slot].channel = channel;
	}

	int nextChild() {
		for(; adjacentVertices.first != adjacentVertices.second; adjacentVertices.first++) {
			const graph_traits<Route::Graph>::vertex_descriptor& child = *(adjacentVertices.first);
			auto incoming = route->getPredecessor(child) == id;
			if(incoming) {
				return child;
			}
		}
		return -1;
	}

	void markChildVisited() {
		adjacentVertices.first++;
	}

	int id;
	int parent;
	int si;
	std::vector<TASCNode*>* nodes;
	TDMASchedule* schedule;
	Route* route;
	std::pair<Route::Graph::adjacency_iterator, Route::Graph::adjacency_iterator> adjacentVertices;
	int maxslots;
};

class TAMCNode : public TASCNode {
public:
	virtual void handleNode(int nextSlot) {
		int v = nextChild();
		if(v != -1) {
			markChildVisited();

			int s = route->getDescendants(v) + 1;
			// avoid first slot
			for(int slot = 1; slot < maxslots; slot++) {
				if(schedule->getNodes()[id].slots[slot].type == TDMASchedule::Type::IDLE) {

					// choose channel
					int c = 11;
					for(; c <= 11+16; c++) {
						if(B.at(slot).find(c) == B.at(slot).end()) {
							break;
						}
					}

					assign(slot,v,TDMASchedule::Type::RX,c);
					blockNeighbors(slot, c, v);

					nodes->at(v)->assign(slot,id,TDMASchedule::Type::TX,c);
					dynamic_cast<TAMCNode*>(nodes->at(v))->blockNeighbors(slot, c, id);

					s--;
					if(s == 0) {
						break;
					}
				}
			}
			auto& node = schedule->getNodes()[id];
			node.printSlots();

			nodes->at(v)->forward(nextSlot);
		}
		else if(id != 0) { // node is leaf or sub tree fully handled
			nodes->at(parent)->backtrack(nextSlot+si,si);
		}
	}

	virtual void blockNeighbors(int slot, int c, int without) {
		auto vs = route->getAdjacentVertices(id);
		for(; vs.first != vs.second; vs.first++) {
			const graph_traits<Route::Graph>::vertex_descriptor& nb = *(vs.first);
			if(nb != without) {
				dynamic_cast<TAMCNode*>(nodes->at(nb))->block(slot,c,true);
			}
		}
	}

	virtual void block(int i, int c, bool forward) {
		B.at(i).insert(c);
		if(forward) {
			dynamic_cast<TAMCNode*>(nodes->at(parent))->block(i,c,false);
		}
	}

	vector<set<int>> B;
};

void TDMAGenerator::createTA(Experiment& experiment, Connections& connections, Route& route, TDMASchedule& schedule, bool multi_channel)
{
	experiment.addIntermediate("slotDuration_us",slotDurationTSCH);

	// Calculate maximum number of slots needed
	int maxslots = 0;
	if(multi_channel) {
		for(int n = 0; n < route.getNodeCount(); n++) {
			maxslots = max(maxslots, requiredSlots(route, n));
		}
	}
	else {
		for(int n = 1; n < route.getNodeCount(); n++) { // without root
			int rxSlotsFromNode = route.getDescendants(n) + 1;
			maxslots += rxSlotsFromNode;
		}
	}

	maxslots += 1; // avoid first slot (for CSMA transmission)

	// Generate empty schedule
	schedule.getNodes().resize(route.getNodeCount());
	for(auto& node : schedule.getNodes()) {
		node.slots.resize(maxslots);
	}

	std::vector<TASCNode*> tanodes;
	for(int n = 0; n < route.getNodeCount(); n++) {
		TASCNode* node;
		if(multi_channel) {
			node = new TAMCNode();
		}
		else {
			node = new TASCNode();
		}
		node->id = n;
		node->parent = route.getPredecessor(n);
		node->nodes = &tanodes;
		node->schedule = &schedule;
		node->adjacentVertices = route.getAdjacentVertices(n);
		node->route = &route;
		node->maxslots = maxslots;
		if(dynamic_cast<TAMCNode*>(node)) {
			dynamic_cast<TAMCNode*>(node)->B.resize(maxslots);
		}
		tanodes.push_back(node);
	}

	tanodes.at(0)->forward(1); // avoid first slot (0)

	for(auto node : tanodes) {
		delete node;
	}
	

	cout << "----------------" << endl;
}

void TDMAGenerator::createOrchestraSBD(Experiment& experiment, Connections& connections, Route& route, TDMASchedule& schedule)
{
	experiment.addIntermediate("slotDuration_us",slotDurationTSCH);

	// Calculate maximum number of slots needed (one extra for common)
	int maxslots = route.getNodeCount()+1;

	// Generate empty schedule
	schedule.getNodes().resize(route.getNodeCount());
	for(auto& node : schedule.getNodes()) {
		node.slots.resize(maxslots);
	}

	cout << "Nodes " << route.getNodeCount() << " " << maxslots << endl;
	for(int n = 0; n < route.getNodeCount(); n++) {
		auto& node = schedule.getNodes()[n];
		int slot = n+1;
		auto parent = route.getPredecessor(n);
		node.slots[slot].type = TDMASchedule::Type::TX;
		node.slots[slot].counterpart = parent;
		schedule.getNodes()[parent].slots[slot].type = TDMASchedule::Type::RX;
		schedule.getNodes()[parent].slots[slot].counterpart = n;
		cout << n << " " << parent << endl;
		node.printSlots();
	}
	cout << "----------------" << endl;
}
