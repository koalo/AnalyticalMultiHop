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

#include "Experiment.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/property_tree/json_parser.hpp>

using namespace std;
using namespace boost::filesystem;
using namespace boost::system;

void Experiment::createFromSetFile(string filename, vector<Experiment>& experiments)
{
	std::ifstream paramsFile(filename);

	if(!paramsFile.is_open()) {
     		cerr << "Could not open parameters file" << endl;
		throw -1;
    	}

	int lineNumber = 0;
	vector<string> header;

	while(!paramsFile.eof()) {
		string line;
		getline(paramsFile, line);

		if(paramsFile.eof()) {
			break;
		}

		if(lineNumber > 0) {
			// create new experiment
			experiments.emplace_back();
		}

		stringstream strstr(line);

		int elementNumber = 0;
		while(strstr) {
			string element;
			strstr >> element;

			if(element == "") {
				continue;
			}

			if(lineNumber == 0) {
				// header
				header.push_back(element);
			}
			else {
				// value
				experiments.back().addParameter(header[elementNumber], element);
			}

			elementNumber++;
		}

		lineNumber++;
	}

	paramsFile.close();
}

void Experiment::addParameter(std::string key, std::string value)
{
	parameters.put(key, value);
}

Topology& Experiment::getTopology()
{
	return topology;
}

Connections& Experiment::getConnections()
{
	return connections;
}

Route& Experiment::getRoute()
{
	return route;
}

RelationSet& Experiment::getRelations()
{
	return relations;
}

TDMASchedule& Experiment::getTDMASchedule()
{
	return tdma_schedule;
}

void Experiment::write(path directory, std::string experiment_file)
{
	// write topology file
	string topology_file = "topology.json";
	topology.write((directory / topology_file).string());
	addIntermediate("topology_file",topology_file);

	// write route file
	string route_file = "route.dot";
	route.write((directory / route_file).string());
	addIntermediate("route_file",route_file);

	// write schedule file
	if(tdma_schedule.isCalculated()) {
		string schedule_file = "schedule.json";
		tdma_schedule.write((directory / schedule_file).string());
		addIntermediate("schedule_file",schedule_file);
	}

	// write experiment file
	boost::property_tree::ptree all;
	all.add_child("parameters", parameters);
	all.add_child("intermediates", intermediates);

	boost::property_tree::json_parser::write_json((directory / experiment_file).string(), all);
}

void Experiment::read(std::string experiment_file)
{
	boost::property_tree::ptree all;
	boost::property_tree::json_parser::read_json(experiment_file, all);

	parameters = all.get_child("parameters");
	intermediates = all.get_child("intermediates");

	// get path relative to experiment file
	boost::filesystem::path f(experiment_file);
	boost::filesystem::path dir = f.parent_path();

	// read route file
	route.read((dir / getIntermediate<string>("route_file")).string());

	// read schedule file
	if(intermediates.count("schedule_file") > 0) {
		auto schedule_file = (dir / getIntermediate<string>("schedule_file")).string();
		if(exists(schedule_file)) {
			tdma_schedule.read(schedule_file);
		}
	}
}

template<>
double Experiment::getParameter<double>(const std::string& key) const {
	std::string s = parameters.get<std::string>(key);
	if(s == "inf") {
		return 1.0/0.0;
	}
	else {
		return parameters.get<double>(key);
	}
}

std::string Experiment::getCSMAResultFileName(std::string experiment_file) {
	// get path relative to experiment file
	boost::filesystem::path f(experiment_file);
	boost::filesystem::path dir = f.parent_path();

	return (dir / "result_csma.json").string();
}

std::string Experiment::getTDMAResultFileName(std::string experiment_file) {
	// get path relative to experiment file
	boost::filesystem::path f(experiment_file);
	boost::filesystem::path dir = f.parent_path();

	return (dir / "result_tdma.json").string();
}
