import networkx as nx
import MDAnalysis
import atom_data as ad
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""
#"
#" Function to parse atom data from atom_data lib
#"
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""
def atomParser(atomData):
	returnDict={}
	for data in atomData["resids"].items():
		for group in data[1].items():
			for atomInfo in group[1].items():
				returnDict[data[0]+":"+atomInfo[0]]={"resname":data[0],"atomName":atomInfo[1]["atomName"], "atomType":atomInfo[1]["atomType"],"charge":atomInfo[1]["charge"],"radius":atomInfo[1]["radius"],"epsilon":atomInfo[1]["epsilon"],"bonds":atomInfo[1]["bonds"],"group":group[0]}
	return returnDict
	
	
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""
#"
#" main function that generate a graph in base of
#" all atoms in a pqr file
#"
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""
def convert(pdb, PSF, ff, parameter_file, output):
	atomData=ad.data(ff, parameter_file)
	graph=nx.Graph()
	u=None
	#getting nodes
	u=MDAnalysis.Universe(pdb)
	sel = None
	try:
		sel=u.select_atoms("protein")
	except:
		sel = u.selectAtoms("protein")
	for atom in sel:
		graph.add_node(str(int(atom.index+1)))
	nodes=graph.nodes()
	##############################################
	##
	## if psf File was provided
	##
	##############################################
	if PSF:
		psf=open(PSF,"r")
		flag=0
		######################################################
		##we will extract conection data and create the graph
		######################################################
		for line in psf:
			#parsing bonds 
			if line=="\n":
				flag=0
			if flag==1:
				splitted=line[:-1].split(" ")
				auxList=[]
				for element in splitted:
					if element!="":
						auxList.append(element)
				i=0
				while i<len(auxList):
					if str(int(auxList[i])) in nodes and str(int(auxList[i+1])) in nodes:
						graph.add_edge(str(int(auxList[i])),str(int(auxList[i+1])))
					i+=2	
			if "!NBOND" in line:
				flag=1
		psf.close()

		#parsing atom data

		psf=open(PSF,"r")
		flag=0
		for line in psf:	
			if line=="\n":
				flag=0	
			if flag==1:
				splitted=line[:-1].split(" ")
				auxList=[]
				for element in splitted:
					if element!="":
						auxList.append(element)
				try:	
					graph.node[auxList[0]]['resname']=auxList[3]
				except:
					pass
				try:					
					graph.node[auxList[0]]['segid']=auxList[1]
				except:
					pass
				try:
					graph.node[auxList[0]]['resnum']=auxList[2]
				except:
					pass
				try:					
					graph.node[auxList[0]]['atomName']=auxList[4]
				except:
					pass
				try:	
					graph.node[auxList[0]]['atomType']=auxList[5]
				except:
					pass
				try:					
					graph.node[auxList[0]]['charge']=auxList[6]
				except:
					pass
				try:
					for item in atomData["resids"][auxList[3]].items():
						for item2 in (item[1].items()):
							if item2[0]==auxList[4]:
								graph.node[auxList[0]]['group']=item[0]
								graph.node[auxList[0]]['epsilon']=item2[1]["epsilon"]
								graph.node[auxList[0]]['radius']=item2[1]["radius"]
				except:
					pass
	
			if "!NATOM" in line:
				flag=1
		psf.close()		
		
	#################################################
	##
	## else, if only pdb file was provided
	##
	#################################################
	else:
		#parsing atom info
		dictAtomData=atomParser(atomData)
		#looping structure
		for residue in sel.residues:
			res=residue.resname
			for atom in residue.atoms:
				if res=="HIS":
					flag=0
					pos=None
					
					try:
						pos=res.atoms.HD1.position
						print pos
						flag+=1
						res="HSD"
					except:
						pass
					
					try:
						pos=res.atoms.HD2.position
						print pos
						flag+=1
						res="HSE"
					except:
						pass
					
					if flag==2: #we found HD1 and HD2 so we have aHSP
						res="HSP"
					elif flag==0:
						res="HSD"						
				try:
					graph.node[str(atom.index+1)]['resname']=  str(atom.resname)
					graph.node[str(atom.index+1)]['segid']=    str(atom.segid)
					graph.node[str(atom.index+1)]['resnum']=   str(residue.resid)
					graph.node[str(atom.index+1)]['atomName']= str(dictAtomData[res+":"+ atom.name]["atomName"])
					graph.node[str(atom.index+1)]['atomType']= str(dictAtomData[res+":"+ atom.name]["atomType"])
					graph.node[str(atom.index+1)]['charge']=   str(dictAtomData[res+":"+ atom.name]["charge"])
					graph.node[str(atom.index+1)]['group']=    str(dictAtomData[res+":"+ atom.name]["group"])
					graph.node[str(atom.index+1)]['epsilon']=  str(dictAtomData[res+":"+ atom.name]["epsilon"])
					graph.node[str(atom.index+1)]['radius']=   str(dictAtomData[res+":"+ atom.name]["radius"])
					graph.node[str(atom.index+1)]['bonds']=""
					aux=""
					for i in range(len(dictAtomData[res+":"+ atom.name]["bonds"])):
						aux+=dictAtomData[res+":"+ atom.name]["bonds"][i]+","
					aux=aux[:-1]
					graph.node[str(atom.index+1)]['bonds']=str(aux)					
				except:
					pass

		#creating peptidic bond
		for res in sel.residues:
			for res2 in sel.residues:
				try:
					if res.atoms.C.segid==res2.atoms.N.segid and res.resnum==res2.resnum-1:
						graph.add_edge(str(res.atoms.C.index+1),str(res2.atoms.N.index+1))
				except:
					pass
		#creating conections
		nodeList=graph.nodes()
		nodeList2=graph.nodes()
		for node in nodeList:
			nodeAtomName=None
			try:
				nodeAtomName=graph.node[node]["atomName"]
				for node2 in nodeList2:
					bondedToNode2=None
					try:
						bondedToNode2=graph.node[node2]["bonds"]
						splitted=bondedToNode2.split(",")
						for atom in splitted:
							if atom==nodeAtomName:
								if graph.node[node]["segid"]==graph.node[node2]["segid"] and graph.node[node]["resname"]==graph.node[node2]["resname"] and graph.node[node]["resnum"]==graph.node[node2]["resnum"]:
									graph.add_edge(node,node2)
					except:
						pass
			except:
				pass
	
	
		

	print "saving topology graph on "+output+"/RIP-MD_Results/topology.gml"
	nx.write_gml(graph,output+"/RIP-MD_Results/topology.gml")
	return graph
