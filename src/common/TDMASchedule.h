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

#ifndef TDMA_SCHEDULE_H
#define TDMA_SCHEDULE_H

#include <string>
#include <vector>

class TDMASchedule {
public:
	enum class Type {IDLE,RX,TX,BLOCKED};

	class Slot {
	public:
		Type type = Type::IDLE;
		int counterpart = 0;
		int channel = 11;
	};

	class Node {
	public:
		std::vector<Slot> slots;
		void printSlots();
	};

	bool isCalculated() {
		return nodes.size() > 0;
	}

	void write(const std::string& filename);
	void read(const std::string& filename);

	std::vector<Node>& getNodes() {
		return nodes;
	}

	unsigned int getTotalTXSlots() {
		return totalTXSlots;
	}

	unsigned int getNodeCount() {
		return nodes.size();
	}

private:
	unsigned int totalTXSlots;
	std::vector<Node> nodes;
};

#endif

