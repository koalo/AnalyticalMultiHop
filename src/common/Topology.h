/*
 * Class for creating topologies of nodes
 *
 * Author:	Florian Meier <florian.meier@koalo.de>
 *		Copyright 2015
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

#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <string>

class TopologyGenerator;

class Topology {
public:
	void write(const std::string& filename);

	std::pair<double,double> getNode(int i);

	struct Frame {
		double xmin, xmax, ymin, ymax;
	};

	Frame getFrame();

	double getDistance(int n1, int n2);

private:
	std::vector<std::pair<double,double> > nodeVector;

	Frame frame;

	friend TopologyGenerator;
};

#endif

