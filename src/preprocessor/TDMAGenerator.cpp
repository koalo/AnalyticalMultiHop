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

#include <algorithm>
#include <random>

using namespace boost;
using namespace std;

static const int aBaseSlotDuration = 60; // symbols
static const int symbolDuration = 16; // us
static const int aNumSuperframeSlots = 16; // slots

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

	void resetChildPointer() {
		adjacentVertices = route->getAdjacentVertices(id);
	}

	int id;
	int parent;
	int si;
	std::vector<TASCNode*>* nodes;
	TDMASchedule* schedule;
	Route* route;
	std::pair<Route::Graph::adjacency_iterator, Route::Graph::adjacency_iterator> adjacentVertices;
	int maxslots;
	int startslot;
	double slotCorrection;
	int handledChildren = 0;
};

class TAMCNode : public TASCNode {
public:
	virtual void searchFreeSlots(int count, int v) {
		for(int slot = startslot; slot < maxslots; slot++) {
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

				count--;
				if(count == 0) {
					break;
				}
			}
		}

		auto& node = schedule->getNodes()[id];
		node.printSlots();

		assert(count == 0);
	}

	virtual void handleNode(int nextSlot) {
		int v = nextChild();
		if(v != -1) {
			markChildVisited();

			int s;
			if(slotCorrection > 0) {
				// fixed size schedule (lSTarget > 0)
				s = 1; // first unconditional slot
			}
			else {
				s = route->getDescendants(v) + 1;
			}

			searchFreeSlots(s,v);
			handledChildren++;

			nodes->at(v)->forward(nextSlot);
		}
		else { // node is leaf or sub tree fully handled
			if(slotCorrection > 0 && handledChildren > 0) {
				// fixed size schedule (lSTarget > 0)
				//int s = floor(slotCorrection*route->getDescendants(id)); // floor, otherwise it is not guaranteed that all nodes find a free slot
				int s = floor((maxslots-startslot)/2);
				if(id == 0) { 
					s = maxslots-startslot;
				}
				s -= handledChildren; // that many children already got a slot
				std::map<int,int> divisor;
				resetChildPointer();
				for(int child = nextChild(); child != -1; markChildVisited(), child = nextChild()) {
					divisor[child] = 0;
				}

				resetChildPointer();
				int elected = nextChild(); // first neighbor (0 would be root)
				cout << "e " << elected << endl;
				while(s > 0) {
					resetChildPointer();
					for(int child = nextChild(); child != -1; markChildVisited(), child = nextChild()) {
						if(divisor[child] * (route->getDescendants(elected)+1) < divisor[elected] * (route->getDescendants(child)+1)) {
							elected = child;
							cout << "e " << elected << endl;
						}
					}

					cout << "elected " << elected << endl;
					searchFreeSlots(1,elected);
					divisor[elected]++;
					s--;
				}
			}

			if(id != 0) { 
				nodes->at(parent)->backtrack(nextSlot+si,si);
			}
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

void TDMAGenerator::createTA(Experiment& experiment, Connections& connections, Route& route, TDMASchedule& schedule, int lSTarget, bool multi_channel, bool tsch, bool cap_reduction, int SO, int MO)
{
	if(tsch) {
		experiment.addIntermediate("slotDuration_us",slotDurationTSCH);
	}
	else {
		int slotDurationDSME = (symbolDuration*aBaseSlotDuration)*(1<<SO);
		experiment.addIntermediate("slotDuration_us",slotDurationDSME);
	}

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

	double slotCorrection = -1;
	if(lSTarget > 0) {
		slotCorrection = lSTarget/(double)maxslots;
		cout << "slotCorrection " << slotCorrection << " " << lSTarget << " " << maxslots << endl;
		maxslots = lSTarget;
	}

	int startslot = 0;
	if(tsch) {
		// avoid first slot (-> 6top)
		maxslots += 1;
		startslot = 1;
	}

	// Generate empty schedule
	schedule.getNodes().resize(route.getNodeCount());
	for(auto& node : schedule.getNodes()) {
		node.slots.resize(maxslots);
	}

	// initialize nodes
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
		node->route = &route;
		node->maxslots = maxslots;
		node->startslot = startslot;
		node->slotCorrection = slotCorrection;
		node->resetChildPointer();

		if(dynamic_cast<TAMCNode*>(node)) {
			dynamic_cast<TAMCNode*>(node)->B.resize(maxslots);
		}
		tanodes.push_back(node);
	}

	tanodes.at(0)->forward(startslot);

	for(auto node : tanodes) {
		delete node;
	}
	

	cout << "----------------" << endl;

	if(!tsch) {
		// Add Beacon and CAP and fill to full superframes
		for(int n = 0; n < route.getNodeCount(); n++) {
			auto& slots = schedule.getNodes()[n].slots;
			int musus = 0;
			if(MO < 0) {
				if(!cap_reduction) {
					musus = ceil(maxslots/7.0);
				}
				else {
					// 8 more slots per superframe,
					// but we have to add the 8 CAP slots
					// at the beginning
					musus = ceil((maxslots+8)/(8.0+7.0));
				}
				// only powers of two are allowed (2**ceil(log2(x)))
				musus = pow(2,ceil(log2(musus)));
			}
			else {
				assert(MO >= SO);
				musus = pow(2,MO-SO);
			}
			
			// Shift already existing slots
			for(int musu = 0; musu < musus; musu++) {
				if(slots.size() > musu*16) {
					auto it = next(slots.begin(),musu*16);
					int additionalSlots = 1; // Beacon
					if(musu == 0 || !cap_reduction) {
						additionalSlots += 8; // CAP
					}
					slots.insert(it,additionalSlots,TDMASchedule::Slot());
				}
			}

			// Fill to full size
			assert(slots.size() <= musus*16);
			slots.resize(musus*16);

			schedule.getNodes()[n].printSlots();
		}
	}
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
