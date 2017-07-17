#!/usr/bin/env python3

import networkx as nx
import numpy as np
import os
import sys
sys.path.insert(0,os.path.dirname(__file__)+'/python')
import analytical_multihop as am

RESULT_DIRECTORY = "results"
WIDTH = 3
HEIGHT = 4
SINK = (0,0)

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

# Iterate over rate
for rate in np.arange(1,2,0.5):
    experiment.intervalUp = rate

    directory = os.path.join(RESULT_DIRECTORY,str(rate))

    # Write experiment files
    experiment.write(directory)

    # Execute the experiment
    result = am.execute(directory, mac='csma')

    # Analyze the result for the node with the highest address
    result = result[max(result.keys())]
    result['rate'] = rate

    print(result)
