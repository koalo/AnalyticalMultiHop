/*
 * Class representing an Experiment
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

#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "Topology.h"
#include "Connections.h"
#include "Route.h"
#include "Relations.h"
#include "TDMASchedule.h"

#include <string>
#include <vector>
#include <map>

#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>

class Experiment {
public:
	static void createFromSetFile(std::string filename, std::vector<Experiment>& experiments);

	template<typename Type>
	Type getParameter(const std::string& key) const {
		return parameters.get<Type>(key);
	}

	template<typename Type>
	Type getIntermediate(const std::string& key) const {
		return intermediates.get<Type>(key);
	}

	template<typename Type>
	void addIntermediate(std::string key, Type value) {
		intermediates.put(key, value);
	}

	Topology& getTopology();
	Connections& getConnections();
	Route& getRoute();
	RelationSet& getRelations();
	TDMASchedule& getTDMASchedule();

	void write(boost::filesystem::path directory, std::string experiment_file);

	void read(std::string experiment_file);

	std::string getResultFileName(std::string experiment_file);

private:
	void addParameter(std::string key, std::string value);

	boost::property_tree::ptree parameters;
	boost::property_tree::ptree intermediates;

	Topology topology;
	Connections connections;
	Route route;
	RelationSet relations;
	TDMASchedule tdma_schedule;

	std::string directory;
};

#endif

