import multiprocessing as mp
import distance
import Angle
import MDAnalysis

def singleParallelDisulfide(parameters):
	i,j = parameters
	toReturn=[]
	try:
		eu_dist=round(distance.euclidean(sel[i].SG.position, sel[j].SG.position),3)
		if eu_dist<=dist:
			v1= sel[i].CB.position-sel[i].SG.position#vector 1 is between a C and the S atoms
			v2= sel[j].CB.position-sel[j].SG.position#same thing
			angle=round(Angle.angle_twoVectors(v1,v2),3)
			if angle>=float(angle1) and angle<=float(angle2):
				node1="CYS:"+sel[i].SG.segid+":"+str(sel[i].SG.resnum)
				node2="CYS:"+sel[j].SG.segid+":"+str(sel[j].SG.resnum)
				edgeName=node1+"-(SS)-"+node2	
				if nodes.has_key(node1) and nodes.has_key(node2):
					toReturn.append(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+str(edgeName)+"\t"+str(round(eu_dist,3))+"\t"+str(round(angle,3))+"\n")				
				
	except:
		pass
	return toReturn
	
def parallelDisulfide(pdbFrameName):
	edges=open(outputFolder+"/RIP-MD_Results/Edges/"+pdbFrameName[0]+".edges","a")
	edges.write("Disulfide Bonds\nSource Node\tTarget Node\tEdge Name\tDistance\tDihedral Angle\n")
	u=MDAnalysis.Universe(pdbFrameName[1])
	cys=u.select_atoms("(resname CYS) and (name CA)")
	cysDict={}
	####################################################################
	##
	## In first place we will create a dictionary, where keys will be
	## each CYS residue. The values will be position
	##
	####################################################################
	for cysAtom in cys: #to create dictionary
		node="CYS:"+cysAtom.segid+":"+str(cysAtom.resid)
		cysDict[node]={"CB":"","SG":"","chain":"","resid":""}
	
	cys=u.select_atoms("(resname CYS) and (name CB SG)")	
	for cysAtom in cys: #to save data in dictionary
		node="CYS:"+cysAtom.segid+":"+str(cysAtom.resid)
		if cysAtom.name=="SG":
			cysDict[node]["SG"]=cysAtom.position
			cysDict[node]["chain"]=cysAtom.segid
			cysDict[node]["resid"]=str(cysAtom.resid)
		if cysAtom.name=="CB":
			cysDict[node]["CB"]=cysAtom.position
			cysDict[node]["chain"]=cysAtom.segid
			cysDict[node]["resid"]=str(cysAtom.resid)
				
	#definition for disulphide is for SG distance of 3A and an angle
	#between 60 and 90 degree
	cysteins=cysDict.items()
	for i in range(len(cysteins)):
		j=i+1
		while j<len(cysteins):
			eu_dist=distance.euclidean(cysteins[i][1]["SG"],cysteins[j][1]["SG"])
			if eu_dist<=dist:
				dihedral_angle=Angle.dihedral(cysteins[i][1]["CB"],cysteins[i][1]["SG"],cysteins[j][1]["CB"],cysteins[j][1]["SG"])
				if dihedral_angle>=float(angle1) and dihedral_angle<=float(angle2):
					node1="CYS:"+cysteins[i][1]["chain"]+":"+cysteins[i][1]["resid"]
					node2="CYS:"+cysteins[j][1]["chain"]+":"+cysteins[j][1]["resid"]
					edgeName=node1+"-(SS)-"+node2	
					if nodes.has_key(node1) and nodes.has_key(node2):
						edges.write(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+str(edgeName)+"\t"+str(round(eu_dist,3))+"\t"+str(round(dihedral_angle,3))+"\n")
			j+=1
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
def parse(pdbDict,Folder, nproc, Nodes, Distance, angle):
	#global variables
	global angle1, angle2, outputFolder, nodes, dist, sel
	nodes=Nodes
	dist=Distance
	outputFolder=Folder
	angles=angle[1:-1].split(",")
	angle1=float(angles[0]) #parsing angle1
	angle2=float(angles[1]) #parsing angle2
	if len(pdbDict)>1:
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing	
		edges=pool.map(parallelDisulfide,(pdbDict.items())) #calling multiprocessing job
	else:
		edges=open(outputFolder+"/RIP-MD_Results/Edges/"+pdbDict.items()[0][0]+".edges","a")
		edges.write("Disulfide Bonds\nSource Node\tTarget Node\tEdge Name\tDistance\tAngle\n")
		u=MDAnalysis.Universe(pdbDict.items()[0][1])
		sel=u.select_atoms("resname CYS").residues
		numberList=[]
		for i in range(len(sel)):
			j=i+1
			while j<len(sel):
				numberList.append([i,j])
				j+=1
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing	
		edgesList=pool.map(singleParallelDisulfide,(numberList)) #calling multiprocessing job
		for edgeList in edgesList: #each processor will return a list with edges, so when we return to the main process we have a list of lists
			for edge in edgeList:
				edges.write(edge)
		edges.write("\n")
		edges.close()		
	return
