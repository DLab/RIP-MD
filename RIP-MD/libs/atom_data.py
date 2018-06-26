def data(ff, parameter):
	#################### extracting data from parameter file ###########################
	pFile=open(parameter, "r")
	flag=0
	parameterDict={}
	for line in pFile:
		if line=="\n":
			flag=0
		if flag==1:
			if line[0]!=" " and line[0]!="!" and line[0]!="\t":
				splitted=line[:-1].split(" ")
				aux=[]
				for item in splitted:
					if item!="":
						aux.append(item)		
				parameterDict[aux[0]]={"epsilon":aux[2],"radius":aux[3]}
		if "V(Lennard-Jones)" in line:
			flag=1
	pFile.close()
	#################### extracting data from topology file ###########################
	residData={"resids":{},"PRES":{}}
	ffFile=open(ff,"r")
	flag=0
	group=-1
	auxDict={}
	actualResid=None
	for line in ffFile:
		if line=="\n":
			flag=0		
		
		#getting resid information
		if flag==1:
			#appending atomdata by group
			if "GROUP" in line:
				group+=1
				#creating group
				residData["resids"][actualResid]["group_"+str(group)]={}
			if "ATOM" ==line[0:4]:
				aux=[]
				data=line[:-1].split(" ")
				for item in data:
					if item!="":
						aux.append(item)
				#appending atom information to corresponding group		
				residData["resids"][actualResid]["group_"+str(group)][aux[1]]={"atomName":aux[1],"atomType":aux[2],"charge":aux[3],"epsilon":None,"radius":None,"bonds":[]}
				try:
					residData["resids"][actualResid]["group_"+str(group)][aux[1]]["epsilon"]=parameterDict[aux[2]]["epsilon"]
					residData["resids"][actualResid]["group_"+str(group)][aux[1]]["radius"]=parameterDict[aux[2]]["radius"]
				except:
					pass
				auxDict[actualResid+":"+aux[1]]="group_"+str(group)
			if "BOND"==line[0:4] or "DOUBLE"==line[0:6]:
				aux=[]
				for item in line[:-1].split(" "):
					if item!="" and item!="BOND" and item!="DOUBLE":
						aux.append(item)
				i=0
				while i<len(aux):
					residData["resids"][actualResid][auxDict[actualResid+":"+aux[i]]][aux[i]]["bonds"].append(aux[i+1])
					i+=2
						
			
		if line[0:4]=="RESI":
			actualResid=(line[:-1].split(" "))[1]
			residData["resids"][actualResid]={}
			flag=1
			group=-1
	ffFile.close()
	#########################################################################################
	##
	## the returned structure has the form
	## {resid:{group:{atom:{atomName:,atomType:,charge:,:epsilon:,radius:,bonds:,}}}}
	##
	#########################################################################################
	return residData
