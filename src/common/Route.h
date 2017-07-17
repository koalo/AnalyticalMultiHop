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

#ifndef ROUTE_H
#define ROUTE_H

#include <boost/graph/adjacency_list.hpp>

class Experiment;
class Connections;
class RouteGenerator;

class Link {
public:
	unsigned int id;
	bool up;
	unsigned int source;
	unsigned int destination;

	unsigned int getAnchor()
	{
		if(up) {
			return source;
		}
		else {
			return destination;
		}
	}
};

class Route {
public:
	typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS,
                       /*Vertex properties=*/boost::property<boost::vertex_distance_t, double>,
                       /*Edge properties=*/boost::property<boost::edge_weight_t, double> >
  		Graph;

	Route();

	Link getLinkByAnchor(int anchor, bool up);
	Link getLinkById(int id);

	int getLinkCount();
	//void setTrafficFactors(int anchor, bool up, double input, double arrival);
	void setBER(int anchor, bool up, double BER);
	//double getInputFactor(int id);
	//double getArrivalFactor(int id);
	double getBER(int id);
	int getDescendants(int node);
	int getDirectChildren(int node);
	int getNodeCount();

	std::pair<Graph::vertex_iterator, Graph::vertex_iterator> getAllVertices();
	std::pair<Graph::adjacency_iterator, Graph::adjacency_iterator> getAdjacentVertices(int node);
	int getPredecessor(int node);

	void write(const std::string& filename);
	void read(const std::string& filename);

private:
	Graph g;

	typedef boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
	typedef boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;

	std::vector<vertex_descriptor> predecessors;
	//std::vector<double> linkInput;
	//std::vector<double> linkArrival;
	//std::vector<double> linkBER;
	std::vector<int> descendants;

	// Since both links for an anchor connect the same nodes
	// and the links are symmetric, we can store a single BER
	// for the anchor.
	std::vector<double> anchorBER;

	int nodes;

	friend RouteGenerator;

	class VertexLabelWriter {
	public:
  		VertexLabelWriter(Route& route)
		: route(route) {
		}

  		template <class Vertex>
  		void operator()(std::ostream& out, const Vertex& v) const {
    			out << "[label=\"" << route.predecessors[v] << " ";
		        out << route.descendants[v] << " ";
		        out << route.anchorBER[v] << "\"]";
  		}

	private:
		Route & route;
	};

	class EdgeLabelWriter {
	public:
  		EdgeLabelWriter(Route& route)
		: route(route) {
		}

  		template <class Edge>
  		void operator()(std::ostream& out, const Edge& e) const {
    			//out << "[label=\"" << e.m_source << " " << e.m_target << "\"]";
    			out << "[]";
  		}

	private:
		Route & route;
	};
};

#endif

