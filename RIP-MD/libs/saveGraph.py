import networkx as nx


########################################################################
##
## this library take one graph and the output folder and save the graph
## in a specific format. available. Edge_List, GEXF, GML, GraphML, JSON
## YAML.
##
########################################################################

def gml(graph,outputfolder):
	nx.write_gml(graph,outputfolder)

def edgelist(graph, outputfolder):
	nx.write_edgelist(graph,outputfolder)
	
def gexf(graph, outputfolder):
	nx.write_gexf(graph,outputfolder)

def graphml(graph, outputfolder):
	nx.write_graphml(graph,outputfolder)

def pajek(graph, outputfolder):
	nx.write_pajek(graph,outputfolder)
	
def saveList(graph, outputfolder):
	nx.write_edgelist(graph,outputfolder)


		
