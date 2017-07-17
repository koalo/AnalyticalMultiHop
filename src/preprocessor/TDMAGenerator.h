/*
 * Class for creating a TDMA schedule
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

#ifndef TDMAGENERATOR_H
#define TDMAGENERATOR_H

#include "Experiment.h"

struct NoFreeSlotException : public std::exception
{
};

class TDMAGenerator {
public:
	static void createTA(Experiment& experiment, Connections& connections, Route& route, TDMASchedule& schedule, bool multi_channel);
	static void createOrchestraSBD(Experiment& experiment, Connections& connections, Route& route, TDMASchedule& schedule);

private:
	static int requiredSlots(Route& route, int node);
};

#endif

