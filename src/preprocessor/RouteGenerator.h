/*
 * Class for creating a static routing tree
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

#ifndef ROUTEGENERATOR_H
#define ROUTEGENERATOR_H

class Experiment;
class Connections;

#include <boost/graph/depth_first_search.hpp>

class RouteGenerator {
public:
	static void create(Experiment& experiment, Connections& connections, Route& route);

	typedef boost::adjacency_list<boost::listS, 
				boost::vecS,
				boost::directedS,
				boost::property<boost::vertex_distance_t, int>
				> TreeGraph;

	//typedef boost::property_map<TreeGraph, boost::vertex_distance_t>::type ChildrenMap;

	class TreeVisitor : public boost::default_dfs_visitor {
	public:
		typedef boost::graph_traits<TreeGraph>::vertex_descriptor vertex_descriptor;

		TreeVisitor(RouteGenerator::TreeGraph& t, std::vector<std::vector<int> >& childrens, std::vector<int>& descendants);
		void finish_vertex(const TreeGraph::vertex_descriptor &s, const TreeGraph &g);


		TreeGraph& t;
		std::vector<std::vector<int> >& childrens;
//		ChildrenMap children_map;
		std::vector<int>& descendants;
	};
};

#endif

