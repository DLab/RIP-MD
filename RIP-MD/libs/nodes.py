#import MDAnalysis
import multiprocessing as mp
import os
import sys
import subprocess
import MDAnalysis
from platform import platform


##################################################
##
## function that write all residues in format
## resname:chain:number
##
##################################################
def write_nodes(outputFolder, pdb):
	output=open(outputFolder+"/RIP-MD_Results/nodes","w")
	output.write("ID\tname\n")
	id_node=1
	nodes={}
	#reading structure
	u = MDAnalysis.Universe(pdb)
	for resid in u.atoms.residues:
		output.write(str(id_node)+"\t"+str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid)+"\n")
		nodes[str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid)]=str(id_node)
		id_node+=1
	return nodes
########################################################
##
## Function to transform  secondary structure from
## eight classes into three classes
##
########################################################
def toEightClass(sec_struct):
	##sect struct from 8 classes to 3 classes
	if sec_struct=="." or sec_struct=="S" or sec_struct=="T" or sec_struct=="C":
		return "C"
	elif sec_struct=="H" or sec_struct=="G" or sec_struct=="I":
		return "H"
	elif sec_struct=="E" or sec_struct=="B":
		return "E"
		
###############################################
##
## function that extract attribute data
##
###############################################
def parallelExtractAttr(parameters):	
	frame, pdbName=parameters	
	dsspName=pdbName[:-3]+"dssp" #generating dssp output name
	#if we are using linux
	flag=0
	if "Linux" in platform():
		command=ripmdPath+"dependencies/linux/dssp/dssp-2.0.4-linux-i386 --i "+pdbName+" -o "+dsspName
		#we will use a subprocess to know if DSSP can not be executed correctly
		proc = subprocess.Popen([command], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
		out, err = proc.communicate()
		if proc.returncode!=0:
			flag=1
	"""////////////////////////////////////////////////////////////////
	// in a close future, implementations for others OS will be
	//implemented, so this if will be extended
	////////////////////////////////////////////////////////////////"""
	out=open(output+"/RIP-MD_Results/Attrs/"+frame+".attr","w") #file to save attr
	out.write("Node ID\tSASA\t8-class sec. struct.\t3-class sec. stuct.\n")		
	u=MDAnalysis.Universe(pdbName)	
	segments=[]
	for segment in u.segments:
		segments.append(segment.segid)
	#if we are not use the correct SO
	if flag!=0:
		reschainDict={}
		dsspDict={}
		for resid in u.atoms.residues:
			reschainDict[str(resid.segment.segid)+":"+str(resid.resid)]=str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid)
			dsspDict[str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid)]=nodes[str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid)]+"\t0\t.\tC"		
		for item in dsspDict.items():
			out.write(item[1]+"\n")
		out.close()	
		return	1
	
	###############################################################
	## now we will parse the DSSP output
	###############################################################	

	File=open(dsspName,"r")
	flag=0 #to know where start to read
	dictDSSP={}
	dssp=[]
	previous=-1 # so allways at the begining it will be true
	actualSegID=0
	for line in File:
		if flag==1:
			if "!" not in line and line !="\n":
				spl=line[:-1].split(" ")
				aux=[]
				for item in spl:
					if item!="":
						aux.append(item)
				if int(aux[1])<previous:
					actualSegID+=1
				previous=int(aux[1])
				aux[2]=segments[actualSegID]
				dssp.append(aux)

		if "#  RESIDUE A" in line:
			flag=1 #in this form, this line will not readed
	reschainDict={}
	dsspDict={}
	for resid in u.atoms.residues:
		reschainDict[str(resid.segment.segid)+":"+str(resid.resid)]=str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid)
		dsspDict[str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid)]=nodes[str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid)]+"\t0\t.\tC"
	try:
		for line in dssp:
			#to see sec struct
			secStruct=line[4]
			if (line[4]!="H") and (line[4]!="B") and (line[4]!="E") and (line[4]!="G") and (line[4]!="I") and (line[4]!="T") and (line[4]!="S"):
				secStruct="."
			cont=0	
			auxCont=0
			for char in line:
				if "," in char and cont==0:
					cont=auxCont
				auxCont+=1
			dsspDict[reschainDict[line[2]+":"+line[1]]]=nodes[reschainDict[line[2]+":"+line[1]]]+"\t"+line[cont-1]+"\t"+secStruct+"\t"+toEightClass(secStruct)
	except:
		reschainDict={}
		dsspDict={}
		for resid in u.atoms.residues:
			reschainDict[str(resid.segment.segid)+":"+str(resid.resid)]=str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid)
			dsspDict[str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid)]=nodes[str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid)]+"\t0\t.\tC"		
		for item in dsspDict.items():
			out.write(item[1]+"\n")
		out.close()	
		return	2	
	#################################
	## saving attribute
	#################################
	for item in dsspDict.items():
		out.write(item[1]+"\n")
	out.close()	
	return 0

	
###############################################
##
## Function that write attributes of nodes
## through time simulation (sec. structure in
## 8 and 3 classes and SASA)
##
###############################################
	
def secStructure(pdbNames, Nodes, outputFolder, nproc, RIPMDPath):
	global nodes, ripmdPath, output, nodeNames, nodeIDs, nodeList
	nodes=Nodes
	ripmdPath=RIPMDPath
	output=outputFolder
	
	nodeFile=open(output+"/RIP-MD_Results/nodes","r")
	nodeList=[]
	for line in nodeFile:
		if "ID" not in line:
			splitted=line[:-1].split("\t")
			nodeList.append(splitted)

	if not os.path.exists(output+"/RIP-MD_Results/Attrs"): #if path does not exist
		os.mkdir(output+"/RIP-MD_Results/Attrs")

	pool=mp.Pool(processes=int(nproc)) #for multiprocessing
	listToWork=[]
	for item in pdbNames.items(): #frame_N, pdb_name
		listToWork.append([item[0],item[1]])
	
	listOfAttr=pool.map(parallelExtractAttr,(listToWork)) #list of attrs will be a list of multiples list with format [ [frame] [ [id of node][sasa][sec struct 8 clases][sec struct 3 classes] ] ]
	flag=0
	for item in listOfAttr:
		if item==1 and flag==0:
			print "Can not detect a DSSP version for your OS. Node atributes will be by default"
			flag=1	
		if item==2 and flag==0:
			print "Error while processing your structure using DSSP. Node atributes will be by default"
			flag=1	
	return
	
	
	
	

