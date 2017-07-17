/*
 * Class for creating topologies of nodes
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

#include "Topology.h"
#include "Experiment.h"

#include <boost/property_tree/json_parser.hpp>

using namespace std;

void Topology::write(const string& filename)
{
	boost::property_tree::ptree tree;
	auto& nodeList = tree.add_child("nodes",boost::property_tree::ptree{});

	for(auto& node : nodeVector) {
		boost::property_tree::ptree item;
		item.put("x",node.first);
		item.put("y",node.second);
		nodeList.push_back(make_pair("",item));
	}

	boost::property_tree::json_parser::write_json(filename, tree);
}

pair<double,double> Topology::getNode(int i)
{
	return nodeVector.at(i);
}

Topology::Frame Topology::getFrame()
{
	return frame;
}

double Topology::getDistance(int n1, int n2)
{
	double x = nodeVector[n1].first - nodeVector[n2].first; 
	double y = nodeVector[n1].second - nodeVector[n2].second; 
	return hypot(x,y);
}

