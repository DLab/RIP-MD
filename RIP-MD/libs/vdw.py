import MDAnalysis
import multiprocessing as mp
import networkx as nx
import math
import distance

##############################################
##
## To calculate coulomb potential for PDB
## in a parallel form (for each snapchot)
##
##############################################

def vdwByFrame(parameters):
	frame, pdb=parameters
	u=MDAnalysis.Universe(pdb)
	edgesFile=open(outputFolder+"/RIP-MD_Results/Edges/"+frame+".edges","a")
	edgesFile.write("vdW contacts\nSource Node\tTarget Node\tEdge Name\tDistance\tpotential (kJ/mole)\n")
	sel=u.select_atoms("protein")
	Range=vdwRange[1:-1].split(",")
	leftRange=float(Range[0])
	rightRange=float(Range[1])
	
	for i in range(len(sel)):
		atom1=str(sel[i].index +1)		
		
		j=i+1
		while j<len(sel):
			atom2=str(sel[j].index +1)
			try:
				euclidean_distance= distance.euclidean(sel[i].position,sel[j].position)
				if euclidean_distance-(float(topology.node[atom1]["radius"])+float(topology.node[atom2]["radius"]))>=leftRange and euclidean_distance-(float(topology.node[atom1]["radius"])+float(topology.node[atom2]["radius"]))<=rightRange:				
					#now we know atoms with a certain distance between atoms radius, so then we will know their number of bonds between them
					#through shortestpath. If they have shortestpath we will asume vdw contact only of shortestpath >=vdwRange.
					#if atoms does not have shortest path (i.e. a double or more chains) we also asume a vdw contact
					cantShortestPath=int(vdwDist)+1
					if int(vdwDist)>0:
						try:
							shortestPath=nx.shortest_path(topology, source=atom1, target=atom2)
							cantShortestPath=len(shortestPath)-1
						except:
							cantShortestPath=int(vdwDist)				
	
					if  cantShortestPath>=int(vdwDist):
						#if these atoms are separed by a number of covalent bonds, we will compute their potential
						#potential is computed as appear in the charmm parameter file V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
						#epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
						#Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
						epsilon=math.sqrt(float(topology.node[atom1]["epsilon"])*float(topology.node[atom2]["epsilon"]))
						RminIJ=(float(topology.node[atom1]["radius"])+float(topology.node[atom2]["radius"]))
						rij=euclidean_distance
						potential=round((4.184*(epsilon*(math.pow((RminIJ/rij),12)-(2*math.pow((RminIJ/rij),6))))),3) #kj/mol
						edgesFile.write(nodes[sel[i].resname+":"+sel[i].segment.segid+":"+str(sel[i].resid)]+"\t"+nodes[sel[j].resname+":"+sel[j].segment.segid+":"+str(sel[j].resid)]+"\t"+sel[i].resname+":"+sel[i].segment.segid+":"+str(sel[i].resid)+":"+sel[i].name+"-(vdW)-"+sel[j].resname+":"+sel[j].segment.segid+":"+str(sel[j].resid)+":"+sel[j].name+"\t"+str(round(euclidean_distance,3))+"\t"+str(potential)+"\n")
			except:
				pass
			j+=1
	edgesFile.write("\n")
	edgesFile.close()
	return

def lj(i):
	Range=vdwRange[1:-1].split(",")
	leftRange=float(Range[0])
	rightRange=float(Range[1])
	listToReturn=[]
	j=i+1
	atom1=str(sel[i].index +1)
	while j<len(atomList):
		atom2=str(sel[j].index +1)
		try:
			euclidean_distance= distance.euclidean(sel[i].position,sel[j].position)
			if euclidean_distance-(float(topology.node[atom1]["radius"])+float(topology.node[atom2]["radius"]))>=leftRange and euclidean_distance-(float(topology.node[atom1]["radius"])+float(topology.node[atom2]["radius"]))<=rightRange:				
				#now we know atoms with a certain distance between atoms radius, so then we will know their number of bonds between them
				#through shortestpath. If they have shortestpath we will asume vdw contact only of shortestpath >=vdwRange.
				#if atoms does not have shortest path (i.e. a double or more chains) we also asume a vdw contact
				cantShortestPath=int(vdwDist)+1
				if int(vdwDist)>0:
					try:
						shortestPath=nx.shortest_path(topology, source=atom1, target=atom2)
						cantShortestPath=len(shortestPath)-1
					except:
						cantShortestPath=int(vdwDist)				

				if  cantShortestPath>=int(vdwDist):
					#if these atoms are separed by a number of covalent bonds, we will compute their potential
					#potential is computed as appear in the charmm parameter file V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
					#epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
					#Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
					epsilon=math.sqrt(float(topology.node[atom1]["epsilon"])*float(topology.node[atom2]["epsilon"]))
					RminIJ=(float(topology.node[atom1]["radius"])+float(topology.node[atom2]["radius"]))
					rij=euclidean_distance
					potential=round((4.184*(epsilon*(math.pow((RminIJ/rij),12)-(2*math.pow((RminIJ/rij),6))))),3) #kj/mole
					#listToReturn.append(str(nodes[str(graph.node[graphNodes[i]]["resname"])+":"+str(graph.node[graphNodes[i]]["chain"])+":"+str(graph.node[graphNodes[i]]["resid"])])+"\t"+str(nodes[str(graph.node[graphNodes[j]]["resname"])+":"+str(graph.node[graphNodes[j]]["chain"])+":"+str(graph.node[graphNodes[j]]["resid"])])+"\t"+str(graph.node[graphNodes[i]]["resname"])+":"+str(graph.node[graphNodes[i]]["chain"])+":"+str(graph.node[graphNodes[i]]["resid"])+":"+str(graph.node[graphNodes[i]]["atomName"])+"-(vdW)-"+str(graph.node[graphNodes[j]]["resname"])+":"+str(graph.node[graphNodes[j]]["chain"])+":"+str(graph.node[graphNodes[j]]["resid"])+":"+str(graph.node[graphNodes[j]]["atomName"])+"\t"+str(euclidean_distance)+"\n")
					listToReturn.append(nodes[sel[i].resname+":"+sel[i].segment.segid+":"+str(sel[i].resid)]+"\t"+nodes[sel[j].resname+":"+sel[j].segment.segid+":"+str(sel[j].resid)]+"\t"+sel[i].resname+":"+sel[i].segment.segid+":"+str(sel[i].resid)+":"+sel[i].name+"-(vdW)-"+sel[j].resname+":"+sel[j].segment.segid+":"+str(sel[j].resid)+":"+sel[j].name+"\t"+str(round(euclidean_distance,3))+"\t"+str(potential)+"\n")
					#print str(round(euclidean_distance,3))+"\t"+potential+"\n"
		except:
			pass
		j+=1
	return listToReturn
	
def parse(args, pdbNames, output, nproc, Nodes, vdw_distance, vdw_range,graphTopo):
	global topology, vdwDist, vdwRange, nodes, outputFolder, pdbFrameName, sel, atomList
	atomList=[]
	topology=graphTopo
	outputFolder=output
	vdwDist=vdw_distance
	nodes=Nodes
	vdwRange=vdw_range
	
	if len(pdbNames)>1:
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		pool.map(vdwByFrame,(pdbNames.items())) #calling in parallel form the vdw function
	else:
		edges=open(output+"/RIP-MD_Results/Edges/"+pdbNames.items()[0][0]+".edges","a")
		edges.write("vdW contacts\nSource Node\tTarget Node\tEdge Name\tDistance\tpotential (kJ/mole)\n")
		u=MDAnalysis.Universe(pdbNames.items()[0][1])
		sel=u.select_atoms("protein")
		for i in range(len(sel)):
			atomList.append(i)
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		edgesList=pool.map(lj,(atomList)) #calling in parallel form the function lj
		for edgeList in edgesList:
			for edge in edgeList:
				edges.write(edge)		
		edges.write("\n")
		edges.close()
					
	return
