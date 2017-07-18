#!/usr/bin/env python3

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
sys.path.insert(0,os.path.dirname(__file__)+'/python')
import analytical_multihop as am

RESULT_DIRECTORY = "results"
WIDTH = 3
HEIGHT = 4
SINK = (0,0)
TESTNODE = (WIDTH-1,HEIGHT-1)

# Generate experiment object
experiment = am.Experiment()

# Generate graph and set positions
G = nx.grid_2d_graph(WIDTH,HEIGHT)
for n in G.node:
    G.node[n]['pos'] = n
experiment.set_graph(G)

# Generate routing tree
possible_predecessors = nx.predecessor(G, source=SINK)
R = nx.DiGraph()
for (k,v) in possible_predecessors.items():
    if len(v) > 0:
        R.add_edge(k,v[0]) # select first
experiment.set_routing(R,SINK)

# Draw graph and routing tree
experiment.draw()

# Generate TDMA schedule
experiment.initialize_schedule(schedule_length = nx.number_of_nodes(G)-1)

for n in G.node:
    if n != SINK:
        # assign the (n-1)th slot to the nth node
        slot_id = experiment.get_index(n)-1
        parent = experiment.get_parent(n)
        print(parent)
        experiment.add_slot(n,slot_id,"TX",parent)
        experiment.add_slot(parent,slot_id,"RX",n)

# Iterate over rate
results = []
for rate in np.arange(0.4,1.2,0.05):
    experiment.intervalUp = 1/rate

    directory = os.path.join(RESULT_DIRECTORY,str(rate))

    # Write experiment files
    experiment.write(directory)

    # Execute the experiment
    for mac in ['CSMA','TDMA']:
        result = am.execute(directory, mac=mac)

        # Analyze the result for a single node
        result = result[experiment.get_index(TESTNODE)]
        result['rate'] = rate
        result['MAC'] = mac

        print(result)
        results.append(result)

# Show and plot the result
df = pd.DataFrame(results)
df.set_index(['rate','MAC'],inplace=True)
df['Rtotal'] *= 100 # get percentage
df = df.unstack()
print(df)
df.plot(y='Rtotal',marker='o')
plt.ylabel("PDR for node "+str(TESTNODE)+" [%]")
plt.xlabel("Rate [1/s]")
plt.show()
