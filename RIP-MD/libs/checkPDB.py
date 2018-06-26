import MDAnalysis

	
def check(Output, pdbList, nproc):
	pdbNames=pdbList
	output=Output
	idList=[]
	repeated=[]
	flag=0
	#we will loop the structure, to discover if exist any problem with the structure
	u = None
	try:
		u = MDAnalysis.Universe(pdbNames["frame_0"])
	except:
		print "Problem while trying to check your structure using MDAnalysis. Maybe you should try to:\n- Use \"protein\" selection\n- Not use --missing_atoms flag\nexiting."
		return 1

	for resid in u.atoms.residues:
		if (str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid)) not in idList:
			idList.append(str(resid.resname)+":"+str(resid.segment.segid)+":"+str(resid.resid))
		else:
			flag=1
			if "Chain "+str(resid.segment.segid)+" Res. Number"+str(resid.resid) not in repeated:
				repeated.append("Chain "+str(resid.segment.segid)+" Res. Number"+str(resid.resid))
	if flag==1:
		print "Some residues are sharing the same number in the same chain ("+str(repeated)[1:-1]+"), exiting"
		return 1
	else:
		return 0
