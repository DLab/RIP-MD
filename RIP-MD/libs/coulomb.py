import multiprocessing as mp
import math
import distance
import MDAnalysis
import networkx as nx

##############################################
##
## To calculate coulomb potential for PDB
## in a parallel form (for each snapchot)
##
##############################################

def potentialByFrame(parameters):
	frame, pdb=parameters
	edges=open(outputFolder+"/RIP-MD_Results/Edges/"+frame+".edges","a")
	edges.write("Coulomb\nSource Node\tTarget Node\tEdge Name\tDistance\tCoulomb Potential (kJ/mole)\n")
	positions={} #dict for position
	u=MDAnalysis.Universe(pdb)
	selection = None
	try:
		selection=u.select_atoms(sel)
	except:
		selection=u.selectAtoms(sel)
	
	for atom in selection:
		positions[str(atom.index+1)] = atom.position	

	try:
		selection=u.select_atoms(sel).residues
	except:
		selection=u.selectAtoms(sel).residues	
		
#	for some in selection:
#		for atom in some.atoms:
#			print atom
				
	for i in range(len(selection)):
		nodeList=[]
		groupList1=[]
		resNode1=selection[i].resname+":"+selection[i].segment.segid+":"+str(selection[i].resid)
	
		#geting node name (that is equal to atomindex +1 cause it begin from 0)
		for atom in selection[i].atoms:
			nodeList.append(str(atom.index+1))
		maxGroup1=0	
		#looking for the number of groups that the residue have
		for node in nodeList:
			try:
				group=int((topology.node[node]["group"].split("_"))[1])
				if group>maxGroup1:
					maxGroup1=group
			except:
				pass
		#if we found more than 1 group we create a list of list to have a local copy of wich atoms form each group
		if maxGroup1!=0: #if the node have "group" attribute
			for k in range(maxGroup1):
				groupList1.append([])
		#looping over atoms to append its to the list of groups
		for node in nodeList:
			try:
				group=int((topology.node[node]["group"].split("_"))[1])
				groupList1[group].append(node) #grouplist in the position group we will puth the node name
			except:
				pass	
		#setting the next residues
		j=i+1
		while j<len(selection):
			resNode2=selection[j].resname+":"+selection[j].segment.segid+":"+str(selection[j].resid)		
			nodeList2=[]
			groupList2=[]		
			#geting node name (that is equal to atomindex +1 cause it begin from 0)
			for atom in selection[j].atoms:
				nodeList2.append(str(atom.index+1))
			maxGroup2=0	
			#looking for the number of groups that the residue have
			for node in nodeList2:
				try:
					group=int((topology.node[node]["group"].split("_"))[1])
					if group>maxGroup2:
						maxGroup2=group
				except:
					pass
			#if we found more than 1 group we create a list of list to have a local copy of wich atoms form each group
			if maxGroup2!=0: #if the node have "group" attribute
				for k in range(maxGroup2):
					groupList2.append([])
			#looping over atoms to append its to the list of groups
			for node in nodeList2:
				try:
					group=int((topology.node[node]["group"].split("_"))[1])
					groupList2[group].append(node) #grouplist in the position group we will puth the node name
				except:
					pass
			
			###########################################
			##
			## now we will compare atoms of differents
			## groups with certain bond distance
			##
			###########################################
			for group in groupList1:
				for group2 in groupList2:
					Vcoulomb=0
					Vdi=0
					Vdd=0
					flag=0				
					atomsFromGroup1=[] #atom names
					distAverage=[]
					#for the i- sumatory
					for node in group:
						atomsFromGroup2=[] #atom names
						atomsFromGroup1.append(topology.node[node]["atomName"])
						# for the j- sumatory
						for node2 in group2:
							atomsFromGroup2.append(topology.node[node2]["atomName"])
							euclidean_distance=distance.euclidean(positions[node],positions[node2])
							if euclidean_distance<=Rrf:	#for the threshold
								distAverage.append(euclidean_distance)
								cantShortestPath=None
								#1-4 potential
								if int(excluded)>0:
									try:
										shortestPath=nx.shortest_path(topology, source=node, target=node2)
										cantShortestPath=len(shortestPath)-1
									except:
										cantShortestPath=int(excluded)				
									#and this shortest path is higher than the coulomb excluded atoms
									if cantShortestPath>=int(excluded):
										#######################################################
										##
										## i,j not excluded, j inside cutoff i
										##
										#######################################################
										Vcoulomb+=(float(topology.node[node]["charge"])*float(topology.node[node2]["charge"]))/euclidean_distance
										flag=1 #to know that we entered here
									if RF:
										Vdd+=((-1*float(topology.node[node]["charge"]))*float(topology.node[node2]["charge"])*Crf*math.pow(euclidean_distance,2))/(2*math.pow(Rrf,3))
										Vdi+=((-1*float(topology.node[node]["charge"]))*float(topology.node[node2]["charge"])*(1-(0.5*Crf)))/Rrf
	
					if flag!=0:
						Vcoulomb*=coulombCte
						Vdd*=coulombCte
						Vdi*=coulombCte
						Vcoulomb=Vcoulomb-(Vdd+Vdi)
						if Vcoulomb>=KbT:
							if nodes.has_key(resNode1) and nodes.has_key(resNode2):
								d=0
								for k in range(len(distAverage)):
									d+=distAverage[k]
								d/=len(distAverage)
								edges.write(str(nodes[resNode1])+"\t"+str(nodes[resNode2])+"\t"+resNode1+":"+str(atomsFromGroup1)+"-(coulomb)-"+resNode2+":"+str(atomsFromGroup2)+"\t"+str(round(d,3))+"\t"+str(Vcoulomb)+"\n")
				j+=1
	edges.write("\n")				
	edges.close()

##############################################
##
## To calculate coulomb potential for PDB
## in a parallel form
##
##############################################
	
def potential(i):
	toReturn=[]
	nodeList=[]
	groupList1=[]
	resNode1=sele[i].resname+":"+sele[i].segment.segid+":"+str(sele[i].resid)

	#geting node name (that is equal to atomindex +1 cause it begin from 0)
	for atom in sele[i].atoms:
		nodeList.append(str(atom.index+1))
	maxGroup1=0	
	#looking for the number of groups that the residue have
	for node in nodeList:
		try:
			group=int((topology.node[node]["group"].split("_"))[1])
			if group>maxGroup1:
				maxGroup1=group
		except:
			pass
	#if we found more than 1 group we create a list of list to have a local copy of wich atoms form each group
	if maxGroup1!=0: #if the node have "group" attribute
		for k in range(maxGroup1):
			groupList1.append([])
	#looping over atoms to append its to the list of groups
	for node in nodeList:
		try:
			group=int((topology.node[node]["group"].split("_"))[1])
			groupList1[group].append(node) #grouplist in the position group we will puth the node name
		except:
			pass	
	#setting the next residues
	j=i+1
	while j<len(sele):
		resNode2=sele[j].resname+":"+sele[j].segment.segid+":"+str(sele[j].resid)		
		nodeList2=[]
		groupList2=[]		
		#geting node name (that is equal to atomindex +1 cause it begin from 0)
		for atom in sele[j].atoms:
			nodeList2.append(str(atom.index+1))
		maxGroup2=0	
		#looking for the number of groups that the residue have
		for node in nodeList2:
			try:
				group=int((topology.node[node]["group"].split("_"))[1])
				if group>maxGroup2:
					maxGroup2=group
			except:
				pass
		#if we found more than 1 group we create a list of list to have a local copy of wich atoms form each group
		if maxGroup2!=0: #if the node have "group" attribute
			for k in range(maxGroup2):
				groupList2.append([])
		#looping over atoms to append its to the list of groups
		for node in nodeList2:
			try:
				group=int((topology.node[node]["group"].split("_"))[1])
				groupList2[group].append(node) #grouplist in the position group we will puth the node name
			except:
				pass

		###########################################
		##
		## now we will compare atoms of differents
		## groups with certain bond distance
		##
		###########################################
		for group in groupList1:
			for group2 in groupList2:
				Vcoulomb=0
				Vdi=0
				Vdd=0
				flag=0				
				atomsFromGroup1=[] #atom names
				distAverage=[]
				#for the i- sumatory
				for node in group:
					atomsFromGroup2=[] #atom names
					atomsFromGroup1.append(topology.node[node]["atomName"])
					# for the j- sumatory
					for node2 in group2:
						atomsFromGroup2.append(topology.node[node2]["atomName"])
						euclidean_distance=distance.euclidean(positionDict[node],positionDict[node2])
						if euclidean_distance<=Rrf:	#for the threshold
							distAverage.append(euclidean_distance)
							cantShortestPath=None
							#1-4 potential
							if int(excluded)>0:
								try:
									shortestPath=nx.shortest_path(topology, source=node, target=node2)
									cantShortestPath=len(shortestPath)-1
								except:
									cantShortestPath=int(excluded)				
								#and this shortest path is higher than the coulomb excluded atoms
								if cantShortestPath>=int(excluded):
									#######################################################
									##
									## i,j not excluded, j inside cutoff i
									##
									#######################################################
									Vcoulomb+=(float(topology.node[node]["charge"])*float(topology.node[node2]["charge"]))/euclidean_distance
									flag=1 #to know that we entered here
								if RF:
									Vdd+=((-1*float(topology.node[node]["charge"]))*float(topology.node[node2]["charge"])*Crf*math.pow(euclidean_distance,2))/(2*math.pow(Rrf,3))
									Vdi+=((-1*float(topology.node[node]["charge"]))*float(topology.node[node2]["charge"])*(1-(0.5*Crf)))/Rrf

				if flag!=0:
					Vcoulomb*=coulombCte
					Vdd*=coulombCte
					Vdi*=coulombCte
					Vcoulomb=Vcoulomb-(Vdd+Vdi)
					if Vcoulomb>=KbT:
						if nodes.has_key(resNode1) and nodes.has_key(resNode2):
							d=0
							for k in range(len(distAverage)):
								d+=distAverage[k]
							d/=len(distAverage)
							toReturn.append(str(nodes[resNode1])+"\t"+str(nodes[resNode2])+"\t"+resNode1+":"+str(atomsFromGroup1)+"-(coulomb)-"+resNode2+":"+str(atomsFromGroup2)+"\t"+str(round(d,3))+"\t"+str(Vcoulomb)+"\n")							
		j+=1

	return toReturn


def parse(pdbNames, output, nproc, Nodes, Topology, Ecs, Erf, T, Krf, Rrf_threshold, coulomb_excluded, reactionField, selection):
	global coulombCte, sele, resNumbers, topology, excluded, positionDict, Rrf,Crf,KbT, nodes, outputFolder, RF, sel
	sel = selection
	RF=reactionField
	outputFolder=output
	nodes=Nodes
	KbT=(8.314487/1000)*float(T)
	Rrf=Rrf_threshold
	positionDict={}
	topology=Topology
	resNumbers=[]
	excluded=coulomb_excluded
	#E0 is equal to 8.8541878 E-12 C**2/JM so we will transform meters to armstrongs and coulomb to elementary charge
	E0=(8.8541878*(math.pow(10,-12)))
	#measures of E0 are C**2/JM with C=6.2415E18 and M=10**-10A so we will have
	#E0=8.8E-12*((6.2415E18)**2/10**-10)e**2/JA
	E0=E0*(math.pow(6.2415*(math.pow(10,18)),2)/math.pow(10,10))
	#coulombCte=1.0/(4*math.pi*E0*Ecs)
	Ke=1389.354566824 #kJ, source http://www.ks.uiuc.edu/Research/namd/2.9/ug/node29.html
	coulombCte=Ke*(1/Ecs)
	Crf=((((2*Ecs)-(2*Erf))*(1+(Krf*Rrf)))-(Erf*math.pow((Krf*Rrf),2)))/(((Ecs+(2*Erf))*(1+(Krf*Rrf)))+(Erf*math.pow((Krf*Rrf),2)))

	#for MD
	if len(pdbNames)>1:
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		edgesList=pool.map(potentialByFrame,(pdbNames.items())) #calling in parallel form the coulomb function
	#for PDB
	else:
		#for a single PDB
		u=MDAnalysis.Universe(pdbNames.items()[0][1])
		sele = None
		try:
			sele=u.select_atoms(sel)
		except:
			sele=u.selectAtoms(sel)
	
		for atom in sele:
			positionDict[str(atom.index+1)] = atom.position

		try:
			sele=u.select_atoms(sel).residues
		except:
			sele=u.selectAtoms(sel).residues
		for i in range(len(sele)):
			resNumbers.append(i)		
		edges=open(output+"/RIP-MD_Results/Edges/"+pdbNames.items()[0][0]+".edges","a")
		edges.write("Coulomb\nSource Node\tTarget Node\tEdge Name\tDistance\tCoulomb Potential (kJ/mole)\n")
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		edgesList=pool.map(potential,(resNumbers)) #calling in parallel form the coulomb function
		for edgeList in edgesList:
			for edge in edgeList:
				edges.write(edge)		
		edges.write("\n")
		edges.close()		

				
	return
	
