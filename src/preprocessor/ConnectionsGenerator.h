/*
 * Class for determining link qualities
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

#ifndef CONNECTIONSGENERATOR_H
#define CONNECTIONSGENERATOR_H

#include <vector>

class Experiment;
class Topology;

class Bin : public std::vector<int> {
public:
	int x, y;
};

class BinBucket : public std::vector<Bin> {
public:
	void init(double xmin, double xmax, double ymin, double ymax, double maxTransmissionRange);
	void push(int id, std::pair<double,double> position);
	int getXBins();
	int getYBins();
	Bin& getBinByCoordinates(int x, int y);

private:
	double binwidth;
	int xbins;
	int ybins;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
};

class ConnectionsGenerator {
public:
	static void create(Experiment& experiment, Topology& topology, Connections& connections);

	static double getReceptionPower(double distance, double Ptx);
	static double getBER(double Prx, double Pn);
	static double getMaximumRange(double Prx, double Pdist);
};

#endif

