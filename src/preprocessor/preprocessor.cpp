/*
 * Program for preprocessing an experiment set
 * to apply the analytical model or a simulation.
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

#include <iostream>
#include <set>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <algorithm>
#include <string>
#include <regex>

#include "Experiment.h"
#include "TopologyGenerator.h"
#include "ConnectionsGenerator.h"
#include "RouteGenerator.h"
#include "TDMAGenerator.h"

using namespace std;
using namespace boost::filesystem;
using namespace boost::system;
namespace po = boost::program_options;

vector<Experiment> experiments;

int main(int argc, char** argv)
{
	/* Read command line arguments */
	po::options_description desc("Usage");
	desc.add_options()
		("help", "produce help message")
		("experiments", po::value<string>()->required(), "tab separated file with the set of experiments")
		("output", po::value<string>()->required(), "output directory")
		("threads", po::value<int>()->default_value(999), "maximum number of threads to use")
		("mac", po::value<string>()->default_value("CSMA"), "TDMA schedule to generate (CSMA for none)")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc) ,vm);

	try {
		po::notify(vm);
	}
	catch(po::error& e)
	{
		cerr << e.what() << endl << endl;
		cerr << desc << endl;
		return 1;
	}

	if(vm.count("help")) {
		cerr << desc << endl;
		return 1;
	}


	/* Read file */
       	Experiment::createFromSetFile(vm["experiments"].as<string>(), experiments);

	
	/* Produce output */
	path output_directory(vm["output"].as<string>());

	for(unsigned int i = 0; i < experiments.size(); i++) {
		auto& experiment = experiments.at(i);

		/* Create topology */
		TopologyGenerator::create(experiment, experiment.getTopology());

		/* Create connections */
		ConnectionsGenerator::create(experiment, experiment.getTopology(), experiment.getConnections());

		/* Create route */
		RouteGenerator::create(experiment, experiment.getConnections(), experiment.getRoute());

		string mac = vm["mac"].as<string>();

		regex dsmeRgx("TAMC_DSME_SO_(.*)_MO_(.*)_CAPRED_(.*)");
		smatch dsmeMatches;

		int lSTarget = -1;
		if(mac != "CSMA") {
			lSTarget = experiment.getParameter<int>("lSTarget");
		}
		transform(mac.begin(), mac.end(), mac.begin(), ::toupper);
		if(mac == "TAMC_TSCH") {
			TDMAGenerator::createTA(experiment, experiment.getConnections(), experiment.getRoute(), experiment.getTDMASchedule(),lSTarget,true,true);
		}
		else if(regex_search(mac,dsmeMatches,dsmeRgx)) {
			int SO = atoi(dsmeMatches[1].str().c_str());
			int MO = atoi(dsmeMatches[2].str().c_str());
			int CAPRED = atoi(dsmeMatches[3].str().c_str());
			cout << "SO " << SO << " MO " << MO << " CAPRED " << CAPRED << endl;
			TDMAGenerator::createTA(experiment, experiment.getConnections(), experiment.getRoute(), experiment.getTDMASchedule(),lSTarget,true,false,CAPRED,SO,MO);
		}
		else if(mac == "TASC_TSCH") {
			TDMAGenerator::createTA(experiment, experiment.getConnections(), experiment.getRoute(), experiment.getTDMASchedule(),lSTarget,false,true);
		}
		else if(mac == "ORCHESTRA") {
			TDMAGenerator::createOrchestraSBD(experiment, experiment.getConnections(), experiment.getRoute(), experiment.getTDMASchedule());
		}
		else if(mac == "CSMA" || mac == "OWN_TSCH") {
			// do nothing
		}
		else {
			cerr << "Invalid mac parameter" << endl;
			cerr << desc << endl;
			return 1;
		}

		/* Produce output */
		// format and create output directory
		boost::format fmt("%03d");
		fmt % i;
		path experiment_directory = output_directory / fmt.str();

		boost::system::error_code err;
		create_directories(experiment_directory, err);
		if(err) {
			cerr << "Could not create output directory!" << endl;
			cerr << experiment_directory << endl;
			return 1;
		}

		// print experiment data
		experiment.write(experiment_directory, vm["mac"].as<string>());
	}
}

