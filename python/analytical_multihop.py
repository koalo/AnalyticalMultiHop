import matplotlib.pyplot as plt
import pygraphviz
import networkx as nx
import os
import json
from networkx.drawing.nx_agraph import write_dot

EXPERIMENT_FILENAME = 'experiment.json'
ROUTE_FILENAME = 'route.dot'

class Experiment: 
    def set_graph(self,G):
        self.G = G

    def set_routing(self,R,sink):
        self.R = R
        self.sink = sink

        # Assign addresses
        self.addresses_for_label = {sink:0}
        self.label_for_address = {0:sink}
        i = 1
        r = nx.reverse(self.R)
        for e in nx.bfs_edges(r,sink):
            if e[1] not in self.addresses_for_label:
                self.addresses_for_label[e[1]] = i
                self.label_for_address[i] = e[1]
                i += 1

    def draw(self):
        nx.draw(self.G,nx.get_node_attributes(self.G,'pos'),edge_color='gray')
        nx.draw(self.R,nx.get_node_attributes(self.G,'pos'),width=4,node_size=2000,labels=self.addresses_for_label)
        plt.show()

    def write(self, result_directory):
        if not os.path.exists(result_directory):
            os.makedirs(result_directory)

        ########################
        # Write experiment file
        experiment = {
            "parameters": {
                "nodes": nx.number_of_nodes(self.G),
                "distance": "25",
                "intervalUp": "0.6666666666666666",
                "intervalDown": "inf",
                "MaxNumberOfBackoffs": "3",
                "MaxNumberOfRetransmissions": "2",
                "InitialBackoffExponent": "3",
                "MaxBackoffExponent": "8",
                "L": "60",
                "handleACKs": "1",
                "betterRetrans": "1",
                "inverse": "0",
                "K": "16",
                "Ptx": "0",
                "Pn": "-80",
                "Pdist": "-80",
                "broadcast": "0"
            },
            "intermediates": {
                "nodesOnOuterCircle": "12",
                "slotDuration_us": "10000",
                "topology_file": "topology.json",
                "route_file": ROUTE_FILENAME,
                "schedule_file": "schedule.json"
            }
        }

        with open(os.path.join(result_directory,EXPERIMENT_FILENAME), 'w') as outfile:
            json.dump(experiment, outfile, indent=4)

        ########################
        # Write route file

        # Generate digraph
        DG = nx.DiGraph(self.G)

        # Clean up
        for n in DG.node:
            DG.node[n] = {}
            #for a in DG.node[n]:
            #    del DG.node[n][a]

        # Assemble label
        DG.node[self.sink]['label'] = "0 "+str(nx.number_of_nodes(self.G)-1)+" 0"
        r = nx.reverse(self.R)
        for n,pred in nx.dfs_predecessors(r,self.sink).items():
            pdesc = len(nx.descendants(r,n))
            DG.node[n]['label'] = str(self.addresses_for_label[pred])+" "+str(pdesc)+" 0"

        # Set addresses
        DG = nx.relabel_nodes(DG, self.addresses_for_label)

        write_dot(DG,os.path.join(result_directory,ROUTE_FILENAME))
