import matplotlib.pyplot as plt
import pygraphviz
import networkx as nx
import os
import json
import subprocess
from networkx.drawing.nx_agraph import write_dot

EXPERIMENT_FILENAME = 'experiment_%s.json'
ROUTE_FILENAME = 'route.dot'
RESULT_FILENAME = "result_%s.json"
SCHEDULE_FILENAME = 'schedule_%s.json'

def execute(result_directory, mac='CSMA'):
    if mac not in ['CSMA','TDMA']:
        raise ValueError("Only csma and tdma allowed for the mac parameter")

    subprocess.call(['./'+mac.lower()+'_model', '--experiment', os.path.join(result_directory,EXPERIMENT_FILENAME%mac)])

    with open(os.path.join(result_directory,RESULT_FILENAME%mac)) as data_file:
            result = json.load(data_file)

    if mac == 'CSMA':
        nodes = len(result)//2+1 # result is given as links!
    else:
        nodes = len(result)

    filtered = {}
    for n in range(1,nodes):
        r = result[str(n)]
        filtered[n] = {}
        filtered[n]['Rtotal'] = float(r['Rtotal'])
        filtered[n]['Dtotal'] = float(r['Dtotal'])

    return filtered

class Experiment: 
    def __init__(self):
        self.intervalUp = 0.6
        self.schedule_length = 0

    def set_graph(self,G):
        self.G = G

    def set_routing(self,R,sink):
        self.R = R
        self.sink = sink

        # Assign indicies
        self.index_for_label = {sink:0}
        self.label_for_index = {0:sink}
        i = 1
        r = nx.reverse(self.R)
        for e in nx.bfs_edges(r,sink):
            if e[1] not in self.index_for_label:
                self.index_for_label[e[1]] = i
                self.label_for_index[i] = e[1]
                i += 1

    def draw(self,block=True):
        nx.draw(self.G,nx.get_node_attributes(self.G,'pos'),edge_color='gray')
        nx.draw(self.R,nx.get_node_attributes(self.G,'pos'),width=4,node_size=2000,labels=self.index_for_label)
        plt.show(block=block)

    def initialize_schedule(self,schedule_length):
        self.schedule_length = schedule_length

        self.schedule = {"nodes":[]}
        for n in self.G.node:
            slots = []
            for i in range(0,schedule_length):
                slots.append({"type": "IDLE",
                              "counterpart": "0",
                              "channel": "0"})
            self.schedule["nodes"].append({"slots":slots})

    def get_index(self,node):
        return self.index_for_label[node]

    def get_node_from_index(self,node_index):
        return self.label_for_index[node_index]

    def get_parent(self,node):
        return self.R.successors(node)[0]

    def add_slot(self,node,slot_id,slot_type,counterpart,channel=0):
        s = self.schedule["nodes"][self.index_for_label[node]]["slots"][slot_id]
        s["type"] = slot_type
        s["counterpart"] = self.index_for_label[counterpart]
        s["channel"] = channel

    def write(self, result_directory, mac):
        if not os.path.exists(result_directory):
            os.makedirs(result_directory)

        ########################
        # Write experiment file
        experiment = {
            "parameters": {
                "nodes": nx.number_of_nodes(self.G),
                "distance": "25",
                "intervalUp": self.intervalUp,
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
                "result_file": RESULT_FILENAME%mac
            }
        }

        if mac != "CSMA":
            experiment["intermediates"]["schedule_file"] = SCHEDULE_FILENAME%mac

        with open(os.path.join(result_directory,EXPERIMENT_FILENAME%mac), 'w') as outfile:
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
            DG.node[n]['label'] = str(self.index_for_label[pred])+" "+str(pdesc)+" 0"

        # Set indicies
        DG = nx.relabel_nodes(DG, self.index_for_label)

        write_dot(DG,os.path.join(result_directory,ROUTE_FILENAME))

        ########################
        # Write schedule file
        if mac != "CSMA":
            assert(self.schedule_length > 0)
            with open(os.path.join(result_directory,SCHEDULE_FILENAME%mac), 'w') as outfile:
                json.dump(self.schedule, outfile, indent=4)
