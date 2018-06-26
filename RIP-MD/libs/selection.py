import MDAnalysis as mda
def save_selection(selection, pdbName, outName):
	try:
		universe = mda.Universe(pdbName)
		#selecting atoms to work
		sel = None
		try:
			sel = universe.selectAtoms(selection) #in v<0.11 or v>0.16 we need to use select_atoms
		except:
			sel = universe.select_atoms(selection) 
		#saving a local copy of selected residues
		writerPDB=mda.Writer(outName)
		writerPDB.write(sel)	
		return outName
	except:
		pass
	return False
