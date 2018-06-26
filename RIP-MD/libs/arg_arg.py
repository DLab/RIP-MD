import multiprocessing as mp
import distance
import MDAnalysis

def singleParallelArgArg(parameters):
	toReturn=[]
	i,j=parameters
	one=ARGs[i].position
	two=ARGs[j].position
	Distance=distance.euclidean(one,two)
	if Distance<=Dist:
		node1=str(ARGs[i].resname)+":"+str(ARGs[i].segment.segid)+":"+str(ARGs[i].resid)
		node2=str(ARGs[j].resname)+":"+str(ARGs[j].segment.segid)+":"+str(ARGs[j].resid)
		edge=node1+"-(Arg-Arg)-"+node2	
		if nodes.has_key(node1) and nodes.has_key(node2):
			toReturn.append(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+str(edge)+"\t"+str(round(Distance,3))+"\n")
	return toReturn


def parallelArgArg(pdbFrameName):
	edges=open(outputFolder+"/RIP-MD_Results/Edges/"+pdbFrameName[0]+".edges","a")
	edges.write("ARG-ARG\nSource Node\tTarget Node\tEdge Name\tDistance (CZ)\n")
	u = MDAnalysis.Universe(pdbFrameName[1])
	#it is based on CZ of ARG residues
	args=u.select_atoms("resname ARG and name CZ") #selecting carbons atoms
	lenArgs=len(args)
	i=0
	j=0
	while i<lenArgs:
		j=i+1
		node1=str(args[i].resname)+":"+str(args[i].segment.segid)+":"+str(args[i].resid)
		while j<lenArgs:
			one=args[i].position
			two=args[j].position
			Distance=distance.euclidean(one,two)
			if Distance<=Dist:
				node2=str(args[j].resname)+":"+str(args[j].segment.segid)+":"+str(args[j].resid)
				edge=node1+"-(Arg-Arg)-"+node2	
				if nodes.has_key(node1) and nodes.has_key(node2):
					edges.write(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+str(edge)+"\t"+str(round(Distance,3))+"\n")
			j+=1
		i+=1
	edges.write("\n")
	edges.close()
	return
		
		
#//////////////////////////////////////////////
#/
#/ main function that will call a parallel
#/ function. this parallel function will
#/ parse and return nodes and attrs for
#/ disulphide bridges
#/
#/////////////////////////////////////////////
def parse(pdbNames,folder, nproc, Nodes, distance):
	global pdbDict, nodes, Dist, outputFolder,ARGs
	
	outputFolder=folder
	nodes=Nodes
	Dist=float(distance)
	if len(pdbNames)>1:
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		pool.map(parallelArgArg,(pdbNames.items()))
	else:
		edges=open(outputFolder+"/RIP-MD_Results/Edges/"+pdbNames.items()[0][0]+".edges","a")
		edges.write("ARG-ARG\nSource Node\tTarget Node\tEdge Name\tDistance (CZ)\n")
		u = MDAnalysis.Universe(pdbNames.items()[0][1])
		ARGs=u.select_atoms("resname ARG and name CZ") #selecting carbons atoms
		parameters=[]
		for i in range(len(ARGs)):
			j=i+1
			while j<len(ARGs):
				parameters.append([i,j])
				j+=1
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		edgesList=pool.map(singleParallelArgArg,(parameters))
		for edgeList in edgesList:
			for edge in edgeList:
				edges.write(edge)		
		edges.write("\n")
		edges.close()					
	return
