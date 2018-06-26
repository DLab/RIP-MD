import multiprocessing as mp
import distance
import MDAnalysis

###############################################
##
## for a single PDB
##
###############################################
def parallelSingleCA(i):
	edgesToReturn=[]
	j=i+1
	one=Calphas[i].position
	node1=str(Calphas[i].resname)+":"+str(Calphas[i].segment.segid)+":"+str(Calphas[i].resid)
	while j<len(Calphas):
		two=Calphas[j].position
		Distance=distance.euclidean(one,two)
		if Distance<=dist:
			node2=str(Calphas[j].resname)+":"+str(Calphas[j].segment.segid)+":"+str(Calphas[j].resid)
			edge=node1+"-(Ca)-"+node2	
			if Nodes.has_key(node1) and Nodes.has_key(node2):
				edgesToReturn.append(str(Nodes[node1])+"\t"+str(Nodes[node2])+"\t"+str(edge)+"\t"+str(round(Distance,3))+"\n")
		j+=1
	return edgesToReturn

###############################################
##
## for MD
##
###############################################

def parallelCA(parameters):
	pdbList=parameters
	edges=open(output+"/RIP-MD_Results/Edges/"+str(pdbList[0])+".edges","a")
	edges.write("Alpha carbon\nSource Node\tTarget Node\tEdge Name\tDistance\n")
	u = MDAnalysis.Universe(pdbList[1])
	calphas=u.select_atoms("name CA") #selecting alpha carbons atoms
	lenCalpha=len(calphas)
	i=0
	j=0
	while i<lenCalpha:
		j=i+1
		node1=str(calphas[i].resname)+":"+str(calphas[i].segment.segid)+":"+str(calphas[i].resid)
		one=calphas[i].position
		while j<lenCalpha:
			two=calphas[j].position
			Distance=distance.euclidean(one,two)
			if Distance<=dist:
				node2=str(calphas[j].resname)+":"+str(calphas[j].segment.segid)+":"+str(calphas[j].resid)
				edge=node1+"-(Ca)-"+node2	
				if Nodes.has_key(node1) and Nodes.has_key(node2):
					edges.write(str(Nodes[node1])+"\t"+str(Nodes[node2])+"\t"+str(edge)+"\t"+str(round(Distance,3))+"\n")
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
#/ Ca distances
#/
#/////////////////////////////////////////////

def parse(pdbList, nproc, Output, nodes, Dist):
	global Nodes, dist,output, Calphas
	Nodes=nodes
	output=Output
	dist=Dist
	listToWork=[]
	if len(pdbList)>1:
		#we will loop over all pdb files
		pdbList=pdbList.items()
		for k in range(len(pdbList)):
			listToWork.append(pdbList[k])
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing			
		edges=pool.map(parallelCA,(listToWork))
	else: #single PDB
		edges=open(output+"/RIP-MD_Results/Edges/"+pdbList.items()[0][0]+".edges","a")
		edges.write("Alpha carbon\nSource Node\tTarget Node\tEdge Name\tDistance\n")
		u = MDAnalysis.Universe(pdbList.items()[0][1])
		calphaN=[]
		Calphas=u.select_atoms("name CA") #selecting alpha carbons atoms
		for i in range(len(Calphas)):
			calphaN.append(i)
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing			
		edgesList=pool.map(parallelSingleCA,(calphaN))
		for edgeList in edgesList: #each processor will return a list with edges, so when we return to the main process we have a list of lists
			for edge in edgeList:
				edges.write(edge)
		edges.write("\n")
		edges.close()	
	return
