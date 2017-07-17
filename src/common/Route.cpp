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
#include "Experiment.h"

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graphviz.hpp>
#include <fstream>

using namespace boost;
using namespace std;

Route::Route() {
}

Link Route::getLinkByAnchor(int anchor, bool up)
{
	Link l;
	l.up = up;

	// example anchor->id for 4 nodes (including 0 as gateway)
	// up:
	// 1 -> 0
	// 2 -> 1
	// 3 -> 2
	//
	// down:
	// 1 -> 3
	// 2 -> 4
	// 3 -> 5

	if(up) {
		l.source = anchor;
		l.destination = predecessors[anchor];
		l.id = anchor-1;
	}
	else {
		l.destination = anchor;
		l.source = predecessors[anchor];
		l.id = anchor+nodes-2;
	}

	return l;
}

Link Route::getLinkById(int id)
{
	Link l;
	l.id = id;
        l.up = id < nodes-1;

	if(l.up) {
		// upstream
		l.source = id+1;
		l.destination = predecessors[l.source];
	}
	else {
		// downstream
		l.destination = id-(nodes-2);
		l.source = predecessors[l.destination];
	}

	return l;
}

int Route::getLinkCount()
{
	return (nodes-1)*2;
}

std::pair<Route::Graph::vertex_iterator, Route::Graph::vertex_iterator> Route::getAllVertices()
{
	return vertices(g);
}

std::pair<Route::Graph::adjacency_iterator, Route::Graph::adjacency_iterator> Route::getAdjacentVertices(int node)
{
	return adjacent_vertices(node, g);
}

int Route::getPredecessor(int node)
{
	return predecessors[node];
}

int Route::getDirectChildren(int node)
{
	return count(predecessors.begin(),predecessors.end(),node);
}

/*
void Route::setTrafficFactors(int anchor, bool up, double input, double arrival)
{
	Link l = getLinkByAnchor(anchor, up);
	linkInput[l.id] = input;
	linkArrival[l.id] = arrival;
}

double Route::getInputFactor(int id)
{
	return linkInput[id];
}

double Route::getArrivalFactor(int id)
{
	return linkArrival[id];
}
*/

double Route::getBER(int id)
{
	//return linkBER[id];
	Link l = getLinkById(id);
	return anchorBER[l.getAnchor()];
}

int Route::getDescendants(int node)
{
	return descendants[node];
}

int Route::getNodeCount() {
	return nodes;
}

void Route::write(const std::string& filename)
{
	ofstream file(filename);
	VertexLabelWriter vlw(*this);
	EdgeLabelWriter elw(*this);
	write_graphviz(file, g, vlw, elw);
}

void Route::read(const std::string& filename)
{
	// Vertex properties
	typedef property < vertex_name_t, std::string,
          property < vertex_color_t, std::string > > vertex_p;

	// Edge properties
	typedef property < edge_weight_t, double > edge_p;

	// adjacency_list-based type
	typedef adjacency_list < vecS, vecS, directedS,
  		vertex_p, edge_p> DotGraph;

	DotGraph dotg;

	dynamic_properties dp;

	property_map<DotGraph, vertex_name_t>::type name = get(vertex_name, dotg);
	dp.property("node_id",name);

	property_map<DotGraph, vertex_color_t>::type color = get(vertex_color, dotg);
	dp.property("label",color);

	// Read graph
	ifstream file(filename);
	read_graphviz(file, dotg, dp, "node_id");

	// Convert graph
	nodes = num_vertices(dotg);
	predecessors.resize(nodes);
	descendants.resize(nodes);
	anchorBER.resize(nodes);

	typedef graph_traits<Graph>::vertex_iterator vertex_iter;
	std::pair<vertex_iter, vertex_iter> vp;
	for (vp = vertices(dotg); vp.first != vp.second; ++vp.first) {
		int id = atoi(get(name, *vp.first).c_str());
		stringstream strstr(get(color, *vp.first));
		strstr >> predecessors[id];
		strstr >> descendants[id];
		strstr >> anchorBER[id];

		add_vertex(*vp.first, g);
	}
	
	typedef graph_traits<DotGraph>::edge_iterator edge_iter;
	std::pair<edge_iter, edge_iter> ep;
	auto iters = edges(dotg);
	edge_iter ei, ei_end;
	std::tie(ei, ei_end) = iters;
	for (; ei != ei_end; ++ei) {
		auto src = source(*ei,dotg);
		int id_src = atoi(get(name, src).c_str());
		auto tgt = target(*ei,dotg);
		int id_tgt = atoi(get(name, tgt).c_str());
		add_edge(id_src, id_tgt, g);
	}
}

