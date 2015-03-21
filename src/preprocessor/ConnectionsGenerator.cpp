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

#include "Connections.h"
#include "ConnectionsGenerator.h"
#include "Experiment.h"
#include "Topology.h"

using namespace std;

void ConnectionsGenerator::create(Experiment& experiment, Topology& topology, Connections& connections)
{
	double Ptx = experiment.getParameter<double>("Ptx");
	double Pn = experiment.getParameter<double>("Pdist");
	double Pdist = experiment.getParameter<double>("Pdist");
	int N = experiment.getParameter<int>("nodes");

	/**
	 * Use bins to speed up the computation:
	 * Do not consider nodes that are too far away.
	 * Reduces the complexity from O(N^2) to approximately O(N*D)
	 * with D being the node density, i.e. the maximum number of
	 * nodes within the range of a node.
	 */
	BinBucket bb;

	double maxTransmissionRange = getMaximumRange(Ptx, Pdist);
	struct Topology::Frame frame = topology.getFrame();
	bb.init(frame.xmin,frame.xmax,frame.ymin,frame.ymax,maxTransmissionRange);

	for(int i = 0; i < N; i++) {
		bb.push(i,topology.getNode(i));
	}

	// reserve some space for the output
	// take 100 as an initial guess for the number of connections per node
	connections.reserve(N*100);

	for(Bin& binA : bb) {
		// iterate over all neighbours, including binA itself (NW,N,NE,W,C,E,SW,S,SE)
		// skip places that are outside the grid
		for(int sx = max(0,binA.x-1); sx <= min(bb.getXBins()-1,binA.x+1); sx++) {
			for(int sy = max(0,binA.y-1); sy <= min(bb.getYBins()-1,binA.y+1); sy++) {
				auto& binB = bb.getBinByCoordinates(sx,sy);

				// iterator over all pair of nodes in binA and binB
				for(auto nA : binA) {
					for(auto nB : binB) {
						if(nA == nB) {
							continue;
						}

						double Prx = getReceptionPower(topology.getDistance(nA,nB), Ptx);

						// allow link only if Prx is high enough to disturb
						if(Prx > Pdist) {
							double BER = getBER(Prx, Pn);
							connections.emplace_back(nA,nB,BER);
						}
					}
				}
			}
		}
	}

	// Only to generate a unambiguous routing tree when comparing to other programs.
	// Can be removed to save time. 
	sort(connections.begin(), connections.end());
}

double ConnectionsGenerator::getReceptionPower(double distance, double Ptx)
{
	double Prx;
	if(distance > 8) {
		Prx = Ptx - (58.5 + 33.0 * log10(distance/8.0));
	}
	else {
		Prx = Ptx - (40.2 + 20.0 * log10(distance));
	}

	return Prx;
}

double ConnectionsGenerator::getMaximumRange(double Ptx, double Pdist)
{
	double receptionPowerAtBreakpoint = getReceptionPower(8, Ptx);
	
	if(receptionPowerAtBreakpoint < Pdist) {
		// reception power at breakpoint is smaller than minimum power needed to disturb
		// -> interference range is smaller than breakpoint distance
		return pow(10,(Ptx-Pdist-40.2)/(20.0));
	}
	else {
		return 8*pow(10,(Ptx-Pdist-58.5)/(33.0));
	}
}

inline int nCk(int n, int k)
{
	int res = 1;
 
	if(k > n-k) {
		k = n - k;
	}
 
	for(int i = 0; i < k; i++) {
		res *= (n - i);
		res /= (i + 1);
	}
 
	return res;
}

double ConnectionsGenerator::getBER(double Prx, double Pn)
{
	double SNR = pow(10,(Prx-Pn)/10.0);

	double s = 0;
	for(int k = 2; k <= 16; k++) {
		s += pow(-1,k)*nCk(16,k)*exp(20.0*SNR*(1.0/k-1));
	}

	return (8.0/15)*(1/16.0)*s;
}

void BinBucket::init(double xmin, double xmax, double ymin, double ymax, double maxTransmissionRange)
{
	this->xmin = xmin;
	this->xmax = xmax;
	this->ymin = ymin;
	this->ymax = ymax;

	// a single bin has minimum width so that a node 
	// on the left edge will not interfer with node on right edge
	binwidth = maxTransmissionRange+1;

	xbins = (xmax-xmin)/binwidth+1;
	ybins = (ymax-ymin)/binwidth+1;

	clear();
	resize(xbins*ybins);

	for(int x = 0; x < xbins; x++) {
		for(int y = 0; y < ybins; y++) {
			Bin& bin = getBinByCoordinates(x,y);
			bin.x = x;
			bin.y = y;
		}
	}
}

Bin& BinBucket::getBinByCoordinates(int x, int y)
{
	return operator[](x*ybins+y);
}

void BinBucket::push(int id, pair<double,double> position)
{
	double bx = (int)((position.first-xmin)/binwidth);
	double by = (int)((position.second-ymin)/binwidth);
	operator[](bx*ybins+by).push_back(id);
}

int BinBucket::getXBins()
{
	return xbins;
}

int BinBucket::getYBins()
{
	return ybins;
}

