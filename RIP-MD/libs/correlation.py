import multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')

from numpy import loadtxt
from os import remove


####################################################
##
## Function to count interactions by node
##
####################################################
def parallelCountInteraction(parameters):
	node, nodeID=parameters
	caFile=open(outputFolder+"/.temp_RIP-MD/"+node+".caByTime","w")	
	hbFile=open(outputFolder+"/.temp_RIP-MD/"+node+".hbByTime","w")	
	saltFile=open(outputFolder+"/.temp_RIP-MD/"+node+".saltByTime","w")	
	ssFile=open(outputFolder+"/.temp_RIP-MD/"+node+".ssByTime","w")	
	cation_piFile=open(outputFolder+"/.temp_RIP-MD/"+node+".cation-piByTime","w")	
	pi_piFile=open(outputFolder+"/.temp_RIP-MD/"+node+".pi-piByTime","w")	
	arg_argFile=open(outputFolder+"/.temp_RIP-MD/"+node+".arg-argByTime","w")	
	coulombFile=open(outputFolder+"/.temp_RIP-MD/"+node+".coulombByTime","w")	
	vdwFile=open(outputFolder+"/.temp_RIP-MD/"+node+".vdwByTime","w")	
	totalFile=open(outputFolder+"/.temp_RIP-MD/"+node+".totalByTime","w")	

	caFile.write(node+"\t"+nodeID+"\n")
	hbFile.write(node+"\t"+nodeID+"\n")
	saltFile.write(node+"\t"+nodeID+"\n")
	ssFile.write(node+"\t"+nodeID+"\n")
	cation_piFile.write(node+"\t"+nodeID+"\n")
	pi_piFile.write(node+"\t"+nodeID+"\n")
	arg_argFile.write(node+"\t"+nodeID+"\n")
	coulombFile.write(node+"\t"+nodeID+"\n")
	vdwFile.write(node+"\t"+nodeID+"\n")
	totalFile.write(node+"\t"+nodeID+"\n")
	for i in range(len(pdbFrameName)):
		edge=open(outputFolder+"/RIP-MD_Results/Edges/frame_"+str(i)+".edges","r")
		countDict={"Ca":0,"HB":0,"salt":0,"SS":0,"cation-pi":0,"pi-pi":0,"Arg-Arg":0,"coulomb":0,"vdW":0,"total":0}
		for line in edge:
			splitted=line.split("\t")
			try:
				if splitted[0]==nodeID or splitted[1]==nodeID:
					if "(Ca)" in line:
						countDict["Ca"]+=1
						countDict["total"]+=1
					elif "(vdW)" in line:
						countDict["vdW"]+=1
						countDict["total"]+=1						
					elif "(salt)" in line:
						countDict["salt"]+=1
						countDict["total"]+=1					
					elif "(Arg-Arg)" in line:
						countDict["Arg-Arg"]+=1
						countDict["total"]+=1						
					elif "(HB)" in line:
						countDict["HB"]+=1
						countDict["total"]+=1						
					elif "(SS)" in line:
						countDict["SS"]+=1
						countDict["total"]+=1						
					elif "(cation-pi)" in line:
						countDict["cation-pi"]+=1
						countDict["total"]+=1						
					elif "(pi-pi)" in line:
						countDict["pi-pi"]+=1
						countDict["total"]+=1						
					elif "(coulomb)" in line:
						countDict["coulomb"]+=1	
						countDict["total"]+=1						
			except: #it is for lines that does not have interactions, like title lines
				pass
													
		caFile.write(str(i)+"\t"+str(countDict["Ca"])+"\n")
		hbFile.write(str(i)+"\t"+str(countDict["HB"])+"\n")
		saltFile.write(str(i)+"\t"+str(countDict["salt"])+"\n")
		ssFile.write(str(i)+"\t"+str(countDict["SS"])+"\n")
		cation_piFile.write(str(i)+"\t"+str(countDict["cation-pi"])+"\n")
		pi_piFile.write(str(i)+"\t"+str(countDict["pi-pi"])+"\n")
		arg_argFile.write(str(i)+"\t"+str(countDict["Arg-Arg"])+"\n")
		coulombFile.write(str(i)+"\t"+str(countDict["coulomb"])+"\n")
		vdwFile.write(str(i)+"\t"+str(countDict["vdW"])+"\n")
		totalFile.write(str(i)+"\t"+str(countDict["total"])+"\n")
		edge.close()
	
	caFile.close()
	hbFile.close()
	saltFile.close()
	ssFile.close()
	cation_piFile.close()
	pi_piFile.close()
	arg_argFile.close()
	coulombFile.close()
	vdwFile.close()
	totalFile.close()

####################################################
##
## Function to do the matrix of interactions where
## a node acts over time
##
####################################################
def matrixNodeTime(pattern):
	
	#counting the number of lines
	X=0
	Y=0
	for node in nodes.items():
		File=open(outputFolder+"/.temp_RIP-MD/"+node[0]+"."+pattern,"r")
		Y=-1
		firstLine=None
		for line in File:
			if Y==-1:
				firstLine=line[:-1]
			if line!="\n":
				Y+=1
		x= int((firstLine.split("\t"))[1])
		if x>X:
			X=x
		File.close()
	#########################################
	## creating a matrix of nodes x time
	#########################################
	matrix=np.zeros((X,Y))
	for node in nodes.items():
		File=open(outputFolder+"/.temp_RIP-MD/"+node[0]+"."+pattern,"r")
		cont=0
		x=None
		for line in File:
			if line!="\n":
				splitted=line[:-1].split("\t")
				if cont==0:
					x=int(splitted[1])-1
					cont+=1
				else:
					y=int(splitted[0])
					value=int(splitted[1])
					matrix[x,y]=value
		File.close()	
	File=open(outputFolder+"/.temp_RIP-MD/"+pattern+".matrix","w")
	for i in range(X):
		for j in range(Y):
			File.write(str(int(matrix[i,j]))+"\t")
		File.write("\n")
	File.close()

####################################################
##
## function to parse the matrix
##
####################################################
def parse(pattern):
	linesFile=[]
	file=open(outputFolder+"/.temp_RIP-MD/"+pattern+".matrix","r")
	for line in file:
		if line!="\n":	
			splitted=line[:-1].split("\t")
			aux=[]
			for num in splitted:
				if num!="":
					aux.append(int(num))
			linesFile.append(aux)
	return linesFile

####################################################
##
## Main function
##
####################################################
def Pearson(pdbNames, output, Nodes, nproc, args):
	global outputFolder, nodes, pdbFrameName
	outputFolder=output
	nodes=Nodes
	pdbFrameName=pdbNames
	
	#we will generate a list of files, where each file have the number of interactions bytime (frame)
	pool=mp.Pool(processes=int(nproc))
	pool.map(parallelCountInteraction,(Nodes.items()))
	
	patterns=[]
	if args.calpha:
		patterns.append("caByTime")
	if args.hbond: ##if hbond option is true
		patterns.append("hbByTime")
	if args.salt:
		patterns.append("saltByTime")
	if args.disulfide:
		patterns.append("ssByTime")
	if args.cation_pi:
		patterns.append("cation-piByTime")	
	if args.pi_pi:
		patterns.append("pi-piByTime")	
	if args.arg_arg:
		patterns.append("arg-argByTime")
	if args.coulomb:
		patterns.append("coulombByTime")
	if args.vdw:
		patterns.append("vdwByTime")
	patterns.append("totalByTime")
			
	pool.map(matrixNodeTime,(patterns))	



	#################################################
	##
	## to compute pearson matrices
	##
	#################################################
	toReturn=[]
	pattern1=0
	while pattern1<len(patterns):
		matrix1=parse(patterns[pattern1])
		pattern2=pattern1
		while pattern2<len(patterns):
			matrix2=parse(patterns[pattern2])
			pearsonMatrix=np.zeros((len(matrix1), len(matrix2)))
			i=0
			while i<len(matrix1):
				j=i
				while j<len(matrix2):
					if i==j:
						pearsonMatrix[i,j]=1
						pearsonMatrix[j,i]=1
					else:
						num=np.corrcoef(matrix1[i],matrix2[j])[0,1]
						if str(num)!="nan":
							pearsonMatrix[i,j]=num
							pearsonMatrix[j,i]=num
						else:
							pearsonMatrix[i,j]=0.0
							pearsonMatrix[j,i]=0.0
						
					j+=1
				i+=1


			aux1=patterns[pattern1].split("By")
			aux2=patterns[pattern2].split("By")
			corr1=open(output+"/RIP-MD_Results/Correlations/"+aux1[0]+"_"+aux2[0],"w")
			corr2=open(output+"/RIP-MD_Results/Correlations/"+aux2[0]+"_"+aux1[0],"w")
			toReturn.append(output+"/RIP-MD_Results/Correlations/"+aux1[0]+"_"+aux2[0])
			toReturn.append(output+"/RIP-MD_Results/Correlations/"+aux2[0]+"_"+aux1[0])
			for i in range(len(matrix1)):
				for j in range(len(matrix2)):
					corr1.write(str(pearsonMatrix[i,j])+"\t")
					corr2.write(str(pearsonMatrix[i,j])+"\t")
				corr1.write("\n")
				corr2.write("\n")
			corr1.close()
			corr2.close()
			pattern2+=1

		pattern1+=1
	
	#some times the diagonal lines have items equal to 0, to assess that this line is 1 we do the following
	file=output+"/RIP-MD_Results/Correlations/"+aux1[0]+"_"+aux2[0]
	matrix=loadtxt(file)
	remove(file)
	file=open(file,"w")
	for i in range(len(matrix[0])):
		for j in range(len(matrix[0])):
			if i==j and matrix[i][j]==0:
				matrix[i][j]=1.0
				print file
			file.write(str(matrix[i][j])+"\t")
		file.write("\n")
	file.close()

	file=output+"/RIP-MD_Results/Correlations/"+aux2[0]+"_"+aux1[0]
	matrix=loadtxt(file)
	remove(file)
	file=open(file,"w")
	for i in range(len(matrix[0])):
		for j in range(len(matrix[0])):
			if i==j and matrix[i][j]==0:
				matrix[i][j]=1.0
				print file
			file.write(str(matrix[i][j])+"\t")
		file.write("\n")
	file.close()	
	return toReturn
