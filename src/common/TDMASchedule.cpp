/*
 * Class for creating TDMA schedules
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

#include "TDMASchedule.h"
#include "Experiment.h"
#include <iostream>
#include <regex>

#include <boost/property_tree/json_parser.hpp>

using namespace std;

void TDMASchedule::write(const string& filename)
{
	boost::property_tree::ptree tree;
	auto& nodeList = tree.add_child("nodes",boost::property_tree::ptree{});

	for(auto& node : nodes) {
		boost::property_tree::ptree nptree;
		auto& slots = nptree.add_child("slots",boost::property_tree::ptree{});
		for(auto& slot : node.slots) {
			boost::property_tree::ptree item;
			switch(slot.type) {
			case Type::IDLE:
				item.put("type","IDLE");
				break;
			case Type::RX:
				item.put("type","RX");
				break;
			case Type::TX:
				item.put("type","TX");
				break;
			case Type::BLOCKED:
				item.put("type","BLOCKED");
				break;
			default:
				assert(false);
			}
			item.put("counterpart",slot.counterpart);
			item.put("channel",slot.channel);
			slots.push_back(make_pair("",item));
		}
		nodeList.push_back(make_pair("",nptree));
	}

	boost::property_tree::json_parser::write_json(filename, tree);
}

void TDMASchedule::read(const std::string& filename)
{
	nodes.clear();

	boost::property_tree::ptree tree;
	boost::property_tree::json_parser::read_json(filename,tree);

	totalTXSlots = 0;

	auto nodes = tree.get_child("nodes");
	for(auto node = nodes.begin(); node != nodes.end(); node++)
	{
		this->nodes.emplace_back();
		auto& nodeObj = this->nodes.back();
		
		auto slots = node->second.get_child("slots");
		for(auto slot = slots.begin(); slot != slots.end(); slot++)
		{
			nodeObj.slots.emplace_back();
			auto& slotObj = nodeObj.slots.back();
			string type = slot->second.get<string>("type");
			if(type == "IDLE") {
				slotObj.type = Type::IDLE;
			}
			else if(type == "BLOCKED") {
				slotObj.type = Type::BLOCKED;
			}
			else {
				if(type == "TX") {
					slotObj.type = Type::TX;
					totalTXSlots++;
				}
				else {
					assert(type == "RX");
					slotObj.type = Type::RX;
				}

				slotObj.counterpart = slot->second.get<int>("counterpart");
				slotObj.channel = slot->second.get<int>("channel");
			}
		}
	}
}

void TDMASchedule::Node::printSlots() {
	for(auto& slot : slots) {
		switch(slot.type) {
		case Type::IDLE:
			std::cout << "I  ";
			break;
		case Type::RX:
			std::cout << "R";
			break;
		case Type::TX:
			std::cout << "T";
			break;
		case Type::BLOCKED:
			std::cout << "B  ";
			break;
		}
		
		if(slot.type == Type::TX || slot.type == Type::RX) {
			std::cout << setfill('0') << setw(2) << slot.counterpart;
		}
		std::cout << " ";
	}
	std::cout << endl;
}

void TDMASchedule::Node::TXfromCommaSeparatedString(const std::string& tx) {
	slots.clear();
	regex re("[\\s,]+");
	sregex_token_iterator it(tx.begin(), tx.end(), re, -1);
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
		slots.push_back(slot);
	}
}
