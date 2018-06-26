import multiprocessing as mp
import numpy as np
import distance
import Angle
import MDAnalysis

############################################
##
## parallel function for a single PDB
##
############################################
def singleParallelCationPi(parameters):
	i,j=parameters
	#i is for the global variable cations
	#j is for piSystems
	toReturn=[]
	###########################################
	#extracting data for the aromatic residue
	###########################################
	aromaticCoords=None
	if piSystems[j].resname=="PHE" or piSystems[j].resname=="TYR":
		try:
			aromaticCoords=[piSystems[j].CG.position,piSystems[j].CD1.position,piSystems[j].CD2.position,piSystems[j].CE1.position,piSystems[j].CE2.position,piSystems[j].CZ.position]
		except:
			pass

	elif piSystems[j].resname=="TRP":
		try:
			aromaticCoords=[piSystems[j].CG.position,piSystems[j].CD1.position,piSystems[j].CD2.position,piSystems[j].NE1.position,piSystems[j].CE2.position,piSystems[j].CE3.position,piSystems[j].CZ2.position,piSystems[j].CZ3.position,piSystems[j].CH2.position]
		except:
			pass
	elif piSystems[j].resname=="HIS" or piSystems[j].resname=="HSD" or piSystems[j].resname=="HSE" or piSystems[j].resname=="HSP":
		try:
			aromaticCoords=[piSystems[j].CG.position,piSystems[j].ND1.position,piSystems[j].CE1.position,piSystems[j].NE2.position,piSystems[j].CD2.position]
		except:
			pass
	if aromaticCoords!=None:
		#getting the centroid
		aux=np.zeros(3)
		for coord in aromaticCoords:
#			if aux==None:
#				aux=coord
#			else:
			aux+=coord
		centroid=aux/len(aromaticCoords)
		###################################
		#extracting charge information
		###################################
		pointOfCharge=np.zeros(3)
		if cations[i].resname=="ARG":
			try:
				pointOfCharge=cations[i].NH1.position
			except:
				pass
		elif cations[i].resname=="LYS":
			try:
				pointOfCharge=cations[i].NZ.position
			except:
				pass
		elif cations[i].resname=="HSP" or cations[i].resname=="HIS":
			try:
				pointOfCharge=cations[i].NE2.postition
			except:
				pass
		if str(pointOfCharge)!=str(np.zeros(3)):
			############################################################
			#now we will compute normal vector for the aromatic residues
			############################################################
			normalVector=np.zeros(3)
			try:
				if piSystems[j].resname=="PHE" or piSystems[j].resname=="TYR":
					v1=piSystems[j].atoms.CG.position-piSystems[j].atoms.CZ.position
					v2=piSystems[j].atoms.CD1.position-piSystems[j].atoms.CE2.position
					v3=piSystems[j].atoms.CD2.position-piSystems[j].atoms.CE1.position
					prod1=np.cross(v1,v2)
					prod2=np.cross(v1,v3)
					prod3=np.cross(v2,v3)
					normalVector=(prod1+prod2+prod3)/3
					
				if piSystems[j].resname=="TRP":
					v1=piSystems[j].atoms.NE1.position-piSystems[j].atoms.CZ3.position
					v2=piSystems[j].atoms.CG.position-piSystems[j].atoms.CH2.position
					v3=piSystems[j].atoms.CD1.position-((piSystems[j].atoms.CZ3.position+piSystems[j].atoms.CH2.position)/2.)
					prod1=np.cross(v1,v2)
					prod2=np.cross(v1,v3)
					prod3=np.cross(v2,v3)
					normalVector=(prod1+prod2+prod3)/3
	
				if piSystems[j].resname=="HIS" or piSystems[j].resname=="HSD" or piSystems[j].resname=="HSE" or piSystems[j].resname=="HSP":
					v1=piSystems[j].atoms.ND1.position-piSystems[j].atoms.NE2.position
					v2=piSystems[j].atoms.CD2.position-piSystems[j].atoms.CE1.position
					v3=piSystems[j].atoms.CG.position-((piSystems[j].atoms.CE1.position+piSystems[j].atoms.NE2.position)/2.)
					prod1=np.cross(v1,v2)
					prod2=np.cross(v1,v3)
					prod3=np.cross(v2,v3)
					normalVector=(prod1+prod2+prod3)/3
			except:
				pass
			if str(normalVector)!=str(np.zeros(3)):
				#computing distance between the charge and the
				dist=distance.euclidean(pointOfCharge,centroid)
				if dist<=Distance: #if distance between cation and center of the ring is less that cutoff
					angle=Angle.angle(pointOfCharge,centroid, normalVector)#calculate the angle between 2 vectors, the vector between centroid of the ring and cation and the normal vector
					if (angle>=a1 and angle<=a2) or (angle<=a3 and angle>=a4):
						node1=None
						node2=None
						try:
							node1=cations[i].resname+":"+cations[i].segment.segid+":"+str(cations[i].resid)
							node2=piSystems[j].resname+":"+piSystems[j].segment.segid+":"+str(piSystems[j].resid)
						except:
							pass
						if nodes.has_key(node1) and nodes.has_key(node2):
							edge_name=node1+"-(cation-pi)-"+node2
							last_less_dist=1000 #a high value to be reemplaced at first loop
							closer_atom=None
							#we will loop over aromatic atoms positions to see what atom is the closer one, with this we will compute the dihedral angle
							for pos in aromaticCoords:
								aux=distance.euclidean(pointOfCharge,pos)
								if aux<last_less_dist:
									last_less_dist=aux
									closer_atom=pos
				#		dihedral angle calculated between CA of the cation, the mass center of cation, closer atom and mass center of ring
						dihedral=None
						try:
							dihedral=Angle.dihedral(cations[i].CA.position,pointOfCharge, closer_atom, centroid)
						except:
							pass
						if dihedral!=None:
							kind=""
							if (dihedral<a5 and dihedral>=a6) or (dihedral<=a7 and dihedral>a8):
								kind="planar"
							elif (dihedral<=a9 and dihedral>=a10)	or (dihedral>=a11 and dihedral<=a12):
								kind="oblique"
							elif (dihedral>a13 and dihedral<a14):
								kind="orthogonal"			
							toReturn.append(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+str(edge_name)+"\t"+str(round(dist,3))+"\t"+str(kind)+"\n")	
				
	return toReturn
	
def parallelCationPi(pdbFrameName):
	edges=open(outputFolder+"/RIP-MD_Results/Edges/"+pdbFrameName[0]+".edges","a")
	edges.write("Cation-Pi\nSource Node\tTarget Node\tEdge Name\tDistance (cation to center of ring)\tOrientation\n")
	u=MDAnalysis.Universe(pdbFrameName[1])
	#selecting residues that act as Cations
	Cations=None
	try:
		Cations=u.select_atoms("resname ARG LYS HSP").residues
	except:
		Cations=u.selectAtoms("resname ARG LYS HSP").residues
	#selecting HIS (resname) with HD1 and HE2 atoms
	try:
		for HIS in u.select_atoms("resname HIS").residues:
			flag=0
			try:
				#we look for the positions of HD1 and HE2, specifics atoms for a protonated HIS,
				#so if this HIS have these atoms we will append the residue to the cation list
				aux=HIS.HD1.position
				flag+=1
			except:
				pass
			try:
				aux=HIS.HE2.position
				flag+=1
			except:
				pass
			if flag==2:
				Cations.append(HIS)
	except:
		for HIS in u.selectAtoms("resname HIS").residues:
			flag=0
			try:
				#we look for the positions of HD1 and HE2, specifics atoms for a protonated HIS,
				#so if this HIS have these atoms we will append the residue to the cation list
				aux=HIS.HD1.position
				flag+=1
			except:
				pass
			try:
				aux=HIS.HE2.position
				flag+=1
			except:
				pass
			if flag==2:
				Cations.append(HIS)		
	#selecting residues that act as a pi system
	PIsystems=None
	try:
		PIsystems=u.select_atoms("resname PHE TRP TYR HIS HSD HSE HSP").residues
	except:
		PIsystems=u.selectAtoms("resname PHE TRP TYR HIS HSD HSE HSP").residues
		
	for i in range(len(Cations)):
		for j in range(len(PIsystems)):
			###########################################
			#extracting data for the aromatic residue
			###########################################
			aromaticCoords=None
			if PIsystems[j].resname=="PHE" or PIsystems[j].resname=="TYR":
				try:
					aromaticCoords=[PIsystems[j].CG.position,PIsystems[j].CD1.position,PIsystems[j].CD2.position,PIsystems[j].CE1.position,PIsystems[j].CE2.position,PIsystems[j].CZ.position]
				except:
					pass
		
			elif PIsystems[j].resname=="TRP":
				try:
					aromaticCoords=[PIsystems[j].CG.position,PIsystems[j].CD1.position,PIsystems[j].CD2.position,PIsystems[j].NE1.position,PIsystems[j].CE2.position,PIsystems[j].CE3.position,PIsystems[j].CZ2.position,PIsystems[j].CZ3.position,PIsystems[j].CH2.position]
				except:
					pass
			elif PIsystems[j].resname=="HIS" or PIsystems[j].resname=="HSD" or PIsystems[j].resname=="HSE" or PIsystems[j].resname=="HSP":
				try:
					aromaticCoords=[PIsystems[j].CG.position,PIsystems[j].ND1.position,PIsystems[j].CE1.position,PIsystems[j].NE2.position,PIsystems[j].CD2.position]
				except:
					pass
			if aromaticCoords!=None:
				#getting the centroid
				aux=np.zeros(3)
				for coord in aromaticCoords:
#					if aux==None:
#						try:
#							print type(coord)
							#aux=coord
							#print aux

					#else:
						#print coord
					aux+=coord
				centroid=aux/len(aromaticCoords)
				###################################
				#extracting charge information
				###################################
				pointOfCharge=np.zeros(3)
				if Cations[i].resname=="ARG":
					try:
						pointOfCharge=Cations[i].NH1.position
					except:
						pass
				elif Cations[i].resname=="LYS":
					try:
						pointOfCharge=Cations[i].NZ.position
					except:
						pass
				elif Cations[i].resname=="HSP" or Cations[i].resname=="HIS":
					try:
						pointOfCharge=Cations[i].NE2.postition
					except:
						pass
				
				if str(pointOfCharge)!=str(np.zeros(3)):
					############################################################
					#now we will compute normal vector for the aromatic residues
					############################################################
					normalVector=np.zeros(3)
					try:
						if PIsystems[j].resname=="PHE" or PIsystems[j].resname=="TYR":
							v1=PIsystems[j].atoms.CG.position-PIsystems[j].atoms.CZ.position
							v2=PIsystems[j].atoms.CD1.position-PIsystems[j].atoms.CE2.position
							v3=PIsystems[j].atoms.CD2.position-PIsystems[j].atoms.CE1.position
							prod1=np.cross(v1,v2)
							prod2=np.cross(v1,v3)
							prod3=np.cross(v2,v3)
							normalVector=(prod1+prod2+prod3)/3
							
						if PIsystems[j].resname=="TRP":
							v1=PIsystems[j].atoms.NE1.position-PIsystems[j].atoms.CZ3.position
							v2=PIsystems[j].atoms.CG.position-PIsystems[j].atoms.CH2.position
							v3=PIsystems[j].atoms.CD1.position-((PIsystems[j].atoms.CZ3.position+PIsystems[j].atoms.CH2.position)/2.)
							prod1=np.cross(v1,v2)
							prod2=np.cross(v1,v3)
							prod3=np.cross(v2,v3)
							normalVector=(prod1+prod2+prod3)/3
			
						if PIsystems[j].resname=="HIS" or PIsystems[j].resname=="HSD" or PIsystems[j].resname=="HSE" or PIsystems[j].resname=="HSP":
							v1=PIsystems[j].atoms.ND1.position-PIsystems[j].atoms.NE2.position
							v2=PIsystems[j].atoms.CD2.position-PIsystems[j].atoms.CE1.position
							v3=PIsystems[j].atoms.CG.position-((PIsystems[j].atoms.CE1.position+PIsystems[j].atoms.NE2.position)/2.)
							prod1=np.cross(v1,v2)
							prod2=np.cross(v1,v3)
							prod3=np.cross(v2,v3)
							normalVector=(prod1+prod2+prod3)/3
					except:
						pass
					if str(normalVector)!=str(np.zeros(3)):
						#computing distance between the charge and the
						dist=distance.euclidean(pointOfCharge,centroid)
						if dist<=Distance: #if distance between cation and center of the ring is less that cutoff
							angle=Angle.angle(pointOfCharge,centroid, normalVector)#calculate the angle between 2 vectors, the vector between centroid of the ring and cation and the normal vector
							if (angle>=a1 and angle<=a2) or (angle<=a3 and angle>=a4):
								node1=None
								node2=None
								try:
									node1=Cations[i].resname+":"+Cations[i].segment.segid+":"+str(Cations[i].resid)
									node2=PIsystems[j].resname+":"+PIsystems[j].segment.segid+":"+str(PIsystems[j].resid)
								except:
									pass
								if nodes.has_key(node1) and nodes.has_key(node2):
									edge_name=node1+"-(cation-pi)-"+node2
									last_less_dist=1000 #a high value to be reemplaced at first loop
									closer_atom=None
									#we will loop over aromatic atoms positions to see what atom is the closer one, with this we will compute the dihedral angle
									for pos in aromaticCoords:
										aux=distance.euclidean(pointOfCharge,pos)
										if aux<last_less_dist:
											last_less_dist=aux
											closer_atom=pos
						#		dihedral angle calculated between CA of the cation, the mass center of cation, closer atom and mass center of ring
								dihedral=None
								try:
									dihedral=Angle.dihedral(Cations[i].CA.position,pointOfCharge, closer_atom, centroid)
								except:
									pass
								if dihedral!=None:
									kind=""
									if (dihedral<a5 and dihedral>=a6) or (dihedral<=a7 and dihedral>a8):
										kind="planar"
									elif (dihedral<=a9 and dihedral>=a10)	or (dihedral>=a11 and dihedral<=a12):
										kind="oblique"
									elif (dihedral>a13 and dihedral<a14):
										kind="orthogonal"			
									edges.write(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+str(edge_name)+"\t"+str(round(dist,3))+"\t"+str(kind)+"\n")	

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
def parse(pdbDict, PSF, folder, nproc, Nodes, dist, angle1, angle2, angle3, angle4):
	global a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14, nodes, Distance, outputFolder, cations, piSystems
	outputFolder=folder
	#####################################
	##
	## first step is parse angles
	##
	#####################################
	#separing angles between cation, centroid and normal vector
	aux=angle1[1:-1].split(",")
	
	Distance=dist
	nodes=Nodes
	a1=float(aux[0])####################     #
	a2=float(aux[1])                   ######## this order is only for my use
	a3=float(aux[3])                   ########
	a4=float(aux[2])####################     #
	#separating angles to define the plane of aromatic ring

	#planar orientation
	aux=angle2[1:-1].split(",")
	a5=float(aux[1]) ####################     #
	a6=float(aux[0])                    ######## this order is only for my use
	a7=float(aux[3])                    ########
	a8=float(aux[2]) ####################     #

	#oblique
	aux=angle3[1:-1].split(",")
	a9=float(aux[1]) ####################     #
	a10=float(aux[0])                   ######## this order is only for my use
	a11=float(aux[2])                   ########
	a12=float(aux[3])####################     #
	
	#orthogonal
	aux=angle4[1:-1].split(",")
	a13=float(aux[0])
	a14=float(aux[1])
	#now we will call a function (with multiprocessing) to get cationpi interactions
	parameters=[]
	#for MD
	if len(pdbDict)>1:
		u=MDAnalysis.Universe(pdbDict.items()[0][1])
		#selecting residues that act as cations
		cations=u.select_atoms("resname ARG LYS HSP").residues
		#selecting HIS (resname) with HD1 and HE2 atoms
		for HIS in u.select_atoms("resname HIS").residues:
			flag=0
			try:
				#we look for the positions of HD1 and HE2, specifics atoms for a protonated HIS,
				#so if this HIS have these atoms we will append the residue to the cation list
				aux=HIS.HD1.position
				flag+=1
			except:
				pass
			try:
				aux=HIS.HE2.position
				flag+=1
			except:
				pass
			if flag==2:
				cations.append(HIS)
		#selecting residues that act as a pi system
		piSystems=u.select_atoms("resname PHE TRP TYR HIS HSD HSE HSP").residues		
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		edges=pool.map(parallelCationPi,(pdbDict.items()))
	#for PDB
	else:
		edges=open(outputFolder+"/RIP-MD_Results/Edges/"+pdbDict.items()[0][0]+".edges","a")
		edges.write("Cation-Pi\nSource Node\tTarget Node\tEdge Name\tDistance (cation to center of ring)\tOrientation\n")
		u=MDAnalysis.Universe(pdbDict.items()[0][1])
		#selecting residues that act as cations
		cations=u.select_atoms("resname ARG LYS HSP").residues
		#selecting HIS (resname) with HD1 and HE2 atoms
		for HIS in u.select_atoms("resname HIS").residues:
			flag=0
			try:
				#we look for the positions of HD1 and HE2, specifics atoms for a protonated HIS,
				#so if this HIS have these atoms we will append the residue to the cation list
				aux=HIS.HD1.position
				flag+=1
			except:
				pass
			try:
				aux=HIS.HE2.position
				flag+=1
			except:
				pass
			if flag==2:
				cations.append(HIS)
		#selecting residues that act as a pi system
		piSystems=u.select_atoms("resname PHE TRP TYR HIS HSD HSE HSP").residues
		parameters=[]
		for i in range(len(cations)):
			for j in range(len(piSystems)):
				parameters.append([i,j])
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		edgesList=pool.map(singleParallelCationPi,(parameters))
		for edgeList in edgesList:
			for edge in edgeList:
				edges.write(edge)
		edges.write("\n")
		edges.close()		
	return
	
