/*
 * Class for determining link qualities
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

#ifndef CONNECTIONS_H
#define CONNECTIONS_H

#include <vector>

class ConnectionsGenerator;

class Connection {
public:
	Connection(int n1, int n2, double BER)
	: n1(n1), n2(n2), BER(BER) {
	}

	bool operator<(const Connection& oth) const {
		return (n1 < oth.n1) || ((n1 == oth.n1) && (n2 < oth.n2));
	}

	int n1, n2;
	double BER;
};

class Connections : private std::vector<Connection> {
public:
	Connection& getConnection(int i);
	int getConnectionCount();

	friend ConnectionsGenerator;
};

#endif

