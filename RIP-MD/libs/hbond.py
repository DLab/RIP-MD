import multiprocessing as mp # to parallelize jobs
import sys

import MDAnalysis
import MDAnalysis.analysis.hbonds as hbonds
import warnings
warnings.filterwarnings("ignore")


def parallelParser(pdbFrameName): #arg: list with universe, start and stop frame, distance and angle to calculate hbond
	sys.stdout = open(outputFolder+"/.temp_RIP-MD/temp_hbond", 'a')				
	sys.stderr = open(outputFolder+"/.temp_RIP-MD/temp_hbond", 'a')
	
	edges_file=open(outputFolder+"/RIP-MD_Results/Edges/"+pdbFrameName[0]+".edges","a") #0 contains the frame_number
	edges_file.write("Hydrogen Bonds\nSource Node\tTarget Node\tEdge Name\tDistance\tAngle\n")
	u = MDAnalysis.Universe(pdbFrameName[1])
	h = hbonds.HydrogenBondAnalysis(u, "protein", "protein", distance=Distance, angle=Angle)
	h.run()
	h.generate_table()
	
	for hbond in h.table:

		node1=str(hbond[5])+":"+dictAtom[hbond[1]]+":"+str(hbond[6])
		node2=str(hbond[8])+":"+dictAtom[hbond[2]]+":"+str(hbond[9])
		if nodes.has_key(node1) and nodes.has_key(node2):
			name=node1+":"+str(hbond[7])+"-(HB)-"+node2+":"+str(hbond[10])
			distance=str(round(float(str(hbond[11])),3))
			angle=str(round(float(str(hbond[12])),3))
			edges_file.write(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+name+"\t"+distance+"\t"+angle+"\n")	
	edges_file.write("\n")
	edges_file.close()

#//////////////////////////////////////////////
#/
#/ main function that will call a parallel
#/ function. this parallel function will
#/ parse and return nodes and attrs for
#/ Hbonds
#/
#/////////////////////////////////////////////
def parse(pdbDict, OutputFolder, nproc, Nodes, RIPMDPath, Dist, Ang):
	global nodes, outputFolder, ripmdPath, Distance, Angle, dictAtom
	nodes=Nodes
	outputFolder=OutputFolder
	ripmdPath=RIPMDPath
	Distance=float(Dist)
	Angle=float(Ang)
	sys.stdout = open(outputFolder+"/.temp_RIP-MD/temp_hbond", 'a')	
	sys.stderr = open(outputFolder+"/.temp_RIP-MD/temp_hbond", 'a')
				
	u = MDAnalysis.Universe(pdbDict["frame_0"])
	dictAtom={}
	chains=[]
	for atom in u.atoms:
		try:
			dictAtom[atom.number]=atom.segid #for versions below 0.16 of mdanalysis
		except:
			dictAtom[atom.id]=atom.segid # for versions higher than 0.16		
	pool=mp.Pool(processes=int(nproc)) #for multiprocessing
	listOfAttr=pool.map(parallelParser,(pdbDict.items()))
	sys.stdout = open(outputFolder+"/output.log", 'a')
	sys.stderr = open(outputFolder+"/output.log", 'a')

	return
