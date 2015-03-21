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

#include "TopologyGenerator.h"
#include "Experiment.h"

#include <boost/property_tree/json_parser.hpp>

using namespace std;

void TopologyGenerator::create(Experiment& experiment, Topology& topology)
{
	int nodes = experiment.getParameter<int>("nodes");
	double mdia = experiment.getParameter<double>("distance");

	double r = 0;
	int circles = 0;

	topology.nodeVector.reserve(nodes);

	// gateway
	nodes--;
	topology.nodeVector.push_back(make_pair(0,0));

	topology.frame.xmin = topology.frame.xmax = topology.frame.ymin = topology.frame.ymax = 0;

	int n = 1;
	while(nodes > 0) {
		r += mdia;

		// for debugging
		// cerr << "On " << circles << " circles: " << nodeVector.size() << " mirrors" << endl;

		double perimeter = 2*M_PI*r;
		int nodesOnPerimeter = perimeter/mdia;

		double angstep = 2*M_PI/nodesOnPerimeter;

		bool lastCircle = false;
		if(nodes <= nodesOnPerimeter) {
			lastCircle = true;
		}

		for(int i = 0; i < nodesOnPerimeter && nodes > 0; i++, nodes--, n++) {
			double ang = i*angstep;
			double x = r*cos(ang);
			double y = r*sin(ang); 
			if(x < topology.frame.xmin) topology.frame.xmin = x;
			if(x > topology.frame.xmax) topology.frame.xmax = x;
			if(y < topology.frame.ymin) topology.frame.ymin = y;
			if(y > topology.frame.ymax) topology.frame.ymax = y;
			topology.nodeVector.push_back(make_pair(x,y));
			if(lastCircle) {
				experiment.addIntermediate("nodesOnOuterCircle", n);
			}
		}

		circles++;
	}
}

