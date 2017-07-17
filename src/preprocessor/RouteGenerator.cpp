/*
 * Class for creating a static routing tree
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

#include "Route.h"
#include "RouteGenerator.h"
#include "Experiment.h"
#include <iostream>

#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;
using namespace std;

RouteGenerator::TreeVisitor::TreeVisitor(RouteGenerator::TreeGraph& t, std::vector<std::vector<int> >& childrens, std::vector<int>& descendants)
: t(t), childrens(childrens), descendants(descendants) {
}

void RouteGenerator::TreeVisitor::finish_vertex(const RouteGenerator::TreeGraph::vertex_descriptor &s, const RouteGenerator::TreeGraph &g) {
	int children = 0;

	// upstream route

	RouteGenerator::TreeGraph::out_edge_iterator edgeItT, edgeEndT;
	boost::tie(edgeItT, edgeEndT) = out_edges( s, t );

	for(; edgeItT != edgeEndT; edgeItT++) {
		vertex_descriptor u = source(*edgeItT, g);
		vertex_descriptor v = target(*edgeItT, g);
		//children += get(children_map, v); // children of child
		children += descendants[v]; // children of child
		children++; // child itself

		if(u == v) { 
			// prevent loops for unconnected nodes          
			cerr << "Unconnected nodes! The results will be inaccurate!" << endl;
			break;
		}

		// downstream route
		childrens[u].push_back(v);

		for(vector<int>::iterator i = childrens[v].begin(); i != childrens[v].end(); i++) {
			int child = *i;
			childrens[u].push_back(child);
		}
	}

	descendants[s] = children;
}

void RouteGenerator::create(Experiment& experiment, Connections& connections, Route& route)
{
	route.nodes = experiment.getParameter<int>("nodes");
	route.predecessors.resize(route.nodes);
	route.descendants.resize(route.nodes);
	route.anchorBER.resize(route.nodes);

	std::map<int, std::map<int, double> > BERs;

	for(int i = 0; i < connections.getConnectionCount(); i++) {
		Connection& link = connections.getConnection(i);

		double minLinkWeight = 1e-3;
		double q = minLinkWeight - log(1 - link.BER);
		
		add_edge(link.n1, link.n2, q, route.g);

		BERs[link.n1][link.n2] = link.BER;
		BERs[link.n2][link.n1] = link.BER;
	}

	route.predecessors.resize(route.nodes);

	graph_traits<Route::Graph>::vertex_descriptor start = vertex(0, route.g);

	dijkstra_shortest_paths(route.g, start,
                          distance_map(get(vertex_distance, route.g)).predecessor_map(boost::make_iterator_property_map(route.predecessors.begin(), get(boost::vertex_index, route.g))));

	/* Set BERs */
	for(int anchor = 1; anchor < route.nodes; anchor++) {
		Link l = route.getLinkByAnchor(anchor, true);
		double BER = BERs[l.source][l.destination];

		// we assume symmetric links
		assert(abs(BER - BERs[l.destination][l.source]) < 10e-15);

		route.anchorBER[l.getAnchor()] = BER;
	}

	/* Walk along tree */
	vector<vector<int> > childrens;
	vector<pair<int,int> > edge_array;
	childrens.resize(route.nodes);
	for(int anchor = 1; anchor < route.nodes; anchor++) {
		Link l = route.getLinkByAnchor(anchor, true);
  		edge_array.push_back(make_pair(l.destination,l.source));
	}

	RouteGenerator::TreeGraph t(edge_array.begin(), edge_array.end(), edge_array.size());

	TreeVisitor vis(t, childrens, route.descendants);
	depth_first_search(t, visitor(vis));
}
