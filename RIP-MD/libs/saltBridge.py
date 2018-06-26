import multiprocessing as mp # to parallelize jobs
import MDAnalysis
import distance

def parallelSalt(pdbFrameName):
	edges=open(outputFolder+"/RIP-MD_Results/Edges/"+pdbFrameName[0]+".edges","a")
	edges.write("Salt bridges\nSource Node\tTarget Node\tEdge Name\tDistance\n")	
	u=MDAnalysis.Universe(pdbFrameName[1])

	# crude definition of salt bridges as contacts between NH/NZ in
	# ARG/LYS and OE*/OD* in ASP/GLU.
	acidic=u.select_atoms("(resname ASP GLU) and (name OE* OD*)")
	basic=u.select_atoms("(resname ARG LYS) and (name NH* NZ)")
	for acidicAtom in acidic:
		acidic_pos=acidicAtom.position
		for basicAtom in basic:
			basic_pos=basicAtom.position
			Distance=distance.euclidean(acidic_pos,basic_pos)
			if Distance <= maxDist:
				#print acidicAtom, basicAtom
				node1=str(acidicAtom.resname)+":"+str(acidicAtom.segment.segid)+":"+str(acidicAtom.resid)
				node2=str(basicAtom.resname)+":"+str(basicAtom.segment.segid)+":"+str(basicAtom.resid)
				name=node1+":"+acidicAtom.name+"-(salt)-"+node2+":"+basicAtom.name
				if nodes.has_key(node1) and nodes.has_key(node2):
					edges.write(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+str(name)+"\t"+str(round(Distance,3))+"\n")
	edges.write("\n")
	edges.close()
	return

#//////////////////////////////////////////////
#/
#/ main function that will call a parallel
#/ function. this parallel function will
#/ parse and return nodes and attrs for
#/ Salt bridges
#/
#/////////////////////////////////////////////
def parse(pdbDict,folder, nproc, Nodes, max_dist):
	global outputFolder, nodes, maxDist
	outputFolder=folder
	nodes=Nodes
	maxDist=max_dist	
	pool=mp.Pool(processes=int(nproc)) #for multiprocessing
	listOfAttr=pool.map(parallelSalt,(pdbDict.items()))
	return
	
