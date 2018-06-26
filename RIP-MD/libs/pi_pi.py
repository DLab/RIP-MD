import multiprocessing as mp
import numpy as np
import distance
import MDAnalysis
import Angle
import math

def singleParallelPiPi(parameters):
	i,j = parameters
	toReturn=[]
	piCoords1=np.zeros(3)
	piCoords2=np.zeros(3)
	centroid1=np.zeros(3)
	centroid2=np.zeros(3)
	normalVector1=np.zeros(3)
	normalVector2=np.zeros(3)
	#geting coords for the first ring
	if piSystems[i].resname=="PHE" or piSystems[i].resname=="TYR":
		try:
			piCoords1=[piSystems[i].atoms.CG.position,piSystems[i].atoms.CD1.position,piSystems[i].atoms.CD2.position,piSystems[i].atoms.CE1.position,piSystems[i].atoms.CE2.position,piSystems[i].atoms.CZ.position]
			#computing normal vector
			vector1=piCoords1[0]-piCoords1[5]
			vector2=piCoords1[1]-piCoords1[4]
			vector3=piCoords1[2]-piCoords1[3]
			prod1=np.cross(vector1,vector2)
			prod2=np.cross(vector1,vector3)
			prod3=np.cross(vector2,vector3)
			normalVector1=(prod1+prod2+prod3)/3
		except:
			pass
	elif piSystems[i].resname=="TRP":
		try:
			piCoords1=[piSystems[i].atoms.CG.position,piSystems[i].atoms.CD1.position,piSystems[i].atoms.CD2.position,piSystems[i].atoms.NE1.position,piSystems[i].atoms.CE2.position,piSystems[i].atoms.CE3.position,piSystems[i].atoms.CZ2.position,piSystems[i].atoms.CZ3.position,piSystems[i].atoms.CH2.position]
			#computing normal vector
			vector1=piCoords1[3]-piCoords1[7]
			vector2=piCoords1[0]-piCoords1[8]
			vector3=piCoords1[1]-((piCoords1[7]+piCoords1[8])/2)
			prod1=np.cross(vector1,vector2)
			prod2=np.cross(vector1,vector3)
			prod3=np.cross(vector2,vector3)
			normalVector1=(prod1+prod2+prod3)/3
		except:
			pass
	elif piSystems[i].resname=="HIS" or piSystems[i].resname=="HSD" or piSystems[i].resname=="HSE" or piSystems[i].resname=="HSP":
		try:
			piCoords1=[piSystems[i].atoms.CG.position,piSystems[i].atoms.ND1.position,piSystems[i].atoms.CE1.position,piSystems[i].atoms.NE2.position,piSystems[i].atoms.CD2.position]
			#computing normal vector			
			vector1=piCoords1[1]-piCoords1[3]
			vector2=piCoords1[0]-piCoords1[2]
			vector3=piCoords1[0]-((piCoords1[2]+piCoords1[3])/2)
			prod1=np.cross(vector1,vector2)
			prod2=np.cross(vector1,vector3)
			prod3=np.cross(vector2,vector3)
			normalVector1=(prod1+prod2+prod3)/3		
		except:
			pass
	#only if the first ring is complete we will compute this with the second one
	if str(piCoords1)!=str(np.zeros(3)):
		for coord in piCoords1:
			if str(centroid1)==str(np.zeros(3)):
				centroid1=coord
			else:
				centroid1+=coord
		centroid1=centroid1/len(piCoords1)
		
		#geting coords for the second ring
		if piSystems[j].resname=="PHE" or piSystems[j].resname=="TYR":
			try:
				piCoords2=[piSystems[j].atoms.CG.position,piSystems[j].atoms.CD1.position,piSystems[j].atoms.CD2.position,piSystems[j].atoms.CE1.position,piSystems[j].atoms.CE2.position,piSystems[j].atoms.CZ.position]
				#computing normal vector
				vector1=piCoords2[0]-piCoords2[5]
				vector2=piCoords2[1]-piCoords2[4]
				vector3=piCoords2[2]-piCoords2[3]
				prod1=np.cross(vector1,vector2)
				prod2=np.cross(vector1,vector3)
				prod3=np.cross(vector2,vector3)
				normalVector2=(prod1+prod2+prod3)/3
			except:
				pass
		elif piSystems[j].resname=="TRP":
			try:
				piCoords2=[piSystems[j].atoms.CG.position,piSystems[j].atoms.CD1.position,piSystems[j].atoms.CD2.position,piSystems[j].atoms.NE1.position,piSystems[j].atoms.CE2.position,piSystems[j].atoms.CE3.position,piSystems[j].atoms.CZ2.position,piSystems[j].atoms.CZ3.position,piSystems[j].atoms.CH2.position]
				#computing normal vector
				vector1=piCoords2[3]-piCoords2[7]
				vector2=piCoords2[0]-piCoords2[8]
				vector3=piCoords2[1]-((piCoords2[7]+piCoords2[8])/2)
				prod1=np.cross(vector1,vector2)
				prod2=np.cross(vector1,vector3)
				prod3=np.cross(vector2,vector3)
				normalVector2=(prod1+prod2+prod3)/3
			except:
				pass
		elif piSystems[j].resname=="HIS" or piSystems[j].resname=="HSD" or piSystems[j].resname=="HSE" or piSystems[j].resname=="HSP":
			try:
				piCoords2=[piSystems[j].atoms.CG.position,piSystems[j].atoms.ND1.position,piSystems[j].atoms.CE1.position,piSystems[j].atoms.NE2.position,piSystems[j].atoms.CD2.position]
				#computing normal vector			
				vector1=piCoords2[1]-piCoords2[3]
				vector2=piCoords2[0]-piCoords2[2]
				vector3=piCoords2[0]-((piCoords2[2]+piCoords2[3])/2)
				prod1=np.cross(vector1,vector2)
				prod2=np.cross(vector1,vector3)
				prod3=np.cross(vector2,vector3)
				normalVector2=(prod1+prod2+prod3)/3
			except:
				pass
		#computing centroid for the second ring
		if str(piCoords2)!=str(np.zeros(3)):
			for coord in piCoords2:
				if str(centroid2)==str(np.zeros(3)):
					centroid2=coord
				else:
					centroid2+=coord
			centroid2/=len(piCoords2)
			
			#computing distance to know if it can be an interaction
			dist=distance.euclidean(centroid1,centroid2)
			if dist<=Dist:
				node1=piSystems[i].resname+":"+piSystems[i].segment.segid+":"+str(piSystems[i].resid)
				node2=piSystems[j].resname+":"+piSystems[j].segment.segid+":"+str(piSystems[j].resid)
				if nodes.has_key(node1) and nodes.has_key(node2):
					edge_name=str(node1)+"-(pi-pi)-"+str(node2)
					#now we define the dihedral angle and  two distances defined as n and p
					#n represents the distance in A between the origin of the orthogonal system (the ring centroid) and the projection of the mass center of the coupled ring j on the Z-axis of the so defined orthogonal system
					#p is the distance in A between the origin of the orthogonal system and the projection of the mass center of the ring j on the XY plane.
					
					#calculating angle between 2 normal vector we can get the dihedral angle
					dihedral=Angle.angle_twoVectors(normalVector1, normalVector2)
					#to get n and p we will change the plane of j-ring and we will use pitagoras theorem
					traslated_centroid=centroid2*1#*1 is to get a copy of coords and not a copy of the pointer to the coords
					traslated_centroid[2]=centroid1[2]
					n=distance.euclidean(centroid1, traslated_centroid)
					p=math.sqrt((centroid2[2]-traslated_centroid[2])**2)
					
					orientation="undefined"
					#now we will define the spatial position
					if (dihedral>=0 and dihedral<30) or (dihedral>=150 and dihedral<180):
						orientation="Parallel orientation"
					
					if (dihedral>=30 and dihedral<150):
						if p<3.5:
							orientation="T-orientation with the edge to face"
						elif p>=3.5 and n<3:
							orientation="T-orientation with the face to the edge"
						else: #p>=3.5 and n >=3
							orientation="L-orientation"
					toReturn.append(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+str(edge_name)+"\t"+str(dist)+"\t"+str(orientation)+"\n")					
	
	return toReturn

def parallelPiPi(pdbFrameName):
	edges=open(outputFolder+"/RIP-MD_Results/Edges/"+pdbFrameName[0]+".edges","a")
	edges.write("Pi-Pi\nSource Node\tTarget Node\tEdge Name\tDistance (center of rings)\tOrientation\n")
	u=MDAnalysis.Universe(pdbFrameName[1])
	##############################################	
	PIsystems=None
	try:
		PIsystems=u.select_atoms("resname PHE TRP TYR HIS HSD HSE HSP").residues
	except:
		PIsystems=u.selectAtoms("resname PHE TRP TYR HIS HSD HSE HSP").residues
		
	for i in range(len(PIsystems)):
		j=i+1
		while j < len(PIsystems):
			#/////////////////////////////////////////////////////////////////////////////////////////////////////////////
			piCoords1=np.zeros(3)
			piCoords2=np.zeros(3)
			centroid1=np.zeros(3)
			centroid2=np.zeros(3)
			normalVector1=np.zeros(3)
			normalVector2=np.zeros(3)
			
			#geting coords for the first ring
			if PIsystems[i].resname=="PHE" or PIsystems[i].resname=="TYR":
				try:
					piCoords1=[PIsystems[i].atoms.CG.position,PIsystems[i].atoms.CD1.position,PIsystems[i].atoms.CD2.position,PIsystems[i].atoms.CE1.position,PIsystems[i].atoms.CE2.position,PIsystems[i].atoms.CZ.position]
					#computing normal vector
					vector1=piCoords1[0]-piCoords1[5]
					vector2=piCoords1[1]-piCoords1[4]
					vector3=piCoords1[2]-piCoords1[3]
					prod1=np.cross(vector1,vector2)
					prod2=np.cross(vector1,vector3)
					prod3=np.cross(vector2,vector3)
					normalVector1=(prod1+prod2+prod3)/3
				except:
					pass
			elif PIsystems[i].resname=="TRP":
				try:
					piCoords1=[PIsystems[i].atoms.CG.position,PIsystems[i].atoms.CD1.position,PIsystems[i].atoms.CD2.position,PIsystems[i].atoms.NE1.position,PIsystems[i].atoms.CE2.position,PIsystems[i].atoms.CE3.position,PIsystems[i].atoms.CZ2.position,PIsystems[i].atoms.CZ3.position,PIsystems[i].atoms.CH2.position]
					#computing normal vector
					vector1=piCoords1[3]-piCoords1[7]
					vector2=piCoords1[0]-piCoords1[8]
					vector3=piCoords1[1]-((piCoords1[7]+piCoords1[8])/2)
					prod1=np.cross(vector1,vector2)
					prod2=np.cross(vector1,vector3)
					prod3=np.cross(vector2,vector3)
					normalVector1=(prod1+prod2+prod3)/3
				except:
					pass
			elif PIsystems[i].resname=="HIS" or PIsystems[i].resname=="HSD" or PIsystems[i].resname=="HSE" or PIsystems[i].resname=="HSP":
				try:
					piCoords1=[PIsystems[i].atoms.CG.position,PIsystems[i].atoms.ND1.position,PIsystems[i].atoms.CE1.position,PIsystems[i].atoms.NE2.position,PIsystems[i].atoms.CD2.position]
					#computing normal vector			
					vector1=piCoords1[1]-piCoords1[3]
					vector2=piCoords1[0]-piCoords1[2]
					vector3=piCoords1[0]-((piCoords1[2]+piCoords1[3])/2)
					prod1=np.cross(vector1,vector2)
					prod2=np.cross(vector1,vector3)
					prod3=np.cross(vector2,vector3)
					normalVector1=(prod1+prod2+prod3)/3		
				except:
					pass
		
			#only if the first ring is complete we will compute this with the second one
			if str(piCoords1)!=str(np.zeros(3)):
				for coord in piCoords1:
					if str(centroid1)==str(np.zeros(3)):
						centroid1=coord
					else:
						centroid1+=coord
				centroid1=centroid1/len(piCoords1)
				
				#geting coords for the second ring
				if PIsystems[j].resname=="PHE" or PIsystems[j].resname=="TYR":
					try:
						piCoords2=[PIsystems[j].atoms.CG.position,PIsystems[j].atoms.CD1.position,PIsystems[j].atoms.CD2.position,PIsystems[j].atoms.CE1.position,PIsystems[j].atoms.CE2.position,PIsystems[j].atoms.CZ.position]
						#computing normal vector
						vector1=piCoords2[0]-piCoords2[5]
						vector2=piCoords2[1]-piCoords2[4]
						vector3=piCoords2[2]-piCoords2[3]
						prod1=np.cross(vector1,vector2)
						prod2=np.cross(vector1,vector3)
						prod3=np.cross(vector2,vector3)
						normalVector2=(prod1+prod2+prod3)/3
					except:
						pass
				elif PIsystems[j].resname=="TRP":
					try:
						piCoords2=[PIsystems[j].atoms.CG.position,PIsystems[j].atoms.CD1.position,PIsystems[j].atoms.CD2.position,PIsystems[j].atoms.NE1.position,PIsystems[j].atoms.CE2.position,PIsystems[j].atoms.CE3.position,PIsystems[j].atoms.CZ2.position,PIsystems[j].atoms.CZ3.position,PIsystems[j].atoms.CH2.position]
						#computing normal vector
						vector1=piCoords2[3]-piCoords2[7]
						vector2=piCoords2[0]-piCoords2[8]
						vector3=piCoords2[1]-((piCoords2[7]+piCoords2[8])/2)
						prod1=np.cross(vector1,vector2)
						prod2=np.cross(vector1,vector3)
						prod3=np.cross(vector2,vector3)
						normalVector2=(prod1+prod2+prod3)/3
					except:
						pass
				elif PIsystems[j].resname=="HIS" or PIsystems[j].resname=="HSD" or PIsystems[j].resname=="HSE" or PIsystems[j].resname=="HSP":
					try:
						piCoords2=[PIsystems[j].atoms.CG.position,PIsystems[j].atoms.ND1.position,PIsystems[j].atoms.CE1.position,PIsystems[j].atoms.NE2.position,PIsystems[j].atoms.CD2.position]
						#computing normal vector			
						vector1=piCoords2[1]-piCoords2[3]
						vector2=piCoords2[0]-piCoords2[2]
						vector3=piCoords2[0]-((piCoords2[2]+piCoords2[3])/2)
						prod1=np.cross(vector1,vector2)
						prod2=np.cross(vector1,vector3)
						prod3=np.cross(vector2,vector3)
						normalVector2=(prod1+prod2+prod3)/3
					except:
						pass
				#computing centroid for the second ring
				if str(piCoords2)!=str(np.zeros(3)):
					for coord in piCoords2:
						if str(centroid2)==str(np.zeros(3)):
							centroid2=coord
						else:
							centroid2+=coord
					centroid2/=len(piCoords2)
					
					#computing distance to know if it can be an interaction
					dist=distance.euclidean(centroid1,centroid2)
					if dist<=Dist:
						node1=PIsystems[i].resname+":"+PIsystems[i].segment.segid+":"+str(PIsystems[i].resid)
						node2=PIsystems[j].resname+":"+PIsystems[j].segment.segid+":"+str(PIsystems[j].resid)
						if nodes.has_key(node1) and nodes.has_key(node2):
							edge_name=str(node1)+"-(pi-pi)-"+str(node2)
							#now we define the dihedral angle and  two distances defined as n and p
							#n represents the distance in A between the origin of the orthogonal system (the ring centroid) and the projection of the mass center of the coupled ring j on the Z-axis of the so defined orthogonal system
							#p is the distance in A between the origin of the orthogonal system and the projection of the mass center of the ring j on the XY plane.
							
							#calculating angle between 2 normal vector we can get the dihedral angle
							dihedral=Angle.angle_twoVectors(normalVector1, normalVector2)
							#to get n and p we will change the plane of j-ring and we will use pitagoras theorem
							traslated_centroid=centroid2*1#*1 is to get a copy of coords and not a copy of the pointer to the coords
							traslated_centroid[2]=centroid1[2]
							n=distance.euclidean(centroid1, traslated_centroid)
							p=math.sqrt((centroid2[2]-traslated_centroid[2])**2)
							
							orientation="undefined"
							#now we will define the spatial position
							if (dihedral>=0 and dihedral<30) or (dihedral>=150 and dihedral<180):
								orientation="Parallel orientation"
							
							if (dihedral>=30 and dihedral<150):
								if p<3.5:
									orientation="T-orientation with the edge to face"
								elif p>=3.5 and n<3:
									orientation="T-orientation with the face to the edge"
								else: #p>=3.5 and n >=3
									orientation="L-orientation"
							edges.write(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+str(edge_name)+"\t"+str(round(dist,3))+"\t"+str(orientation)+"\n")					
			
			#/////////////////////////////////////////////////////////////////////////////////////////////////////////////				
			j+=1

	###############################################
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
	global pdbDict, nodes, Dist, outputFolder, piSystems
	
	outputFolder=folder
	nodes=Nodes
	Dist=distance
	if len(pdbNames)>1:
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		pool.map(parallelPiPi,(pdbNames.items()))
	else:
		edges=open(outputFolder+"/RIP-MD_Results/Edges/"+pdbNames.items()[0][0]+".edges","a")
		edges.write("Pi-Pi\nSource Node\tTarget Node\tEdge Name\tDistance (center of rings)\tOrientation\n")
		u=MDAnalysis.Universe(pdbNames.items()[0][1])
		piSystems = None
		try:
			piSystems=u.select_atoms("resname PHE TRP TYR HIS HSD HSE HSP").residues
		except:
			piSystems=u.selectAtoms("resname PHE TRP TYR HIS HSD HSE HSP").residues
		parameters=[]
		for i in range(len(piSystems)):
			j=i+1
			while j < len(piSystems):
				parameters.append([i,j])
				j+=1
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		edgesList=pool.map(singleParallelPiPi,(parameters))
		for edgeList in edgesList:
			for edge in edgeList:
				edges.write(edge)
		edges.write("\n")
		edges.close()		
		return
