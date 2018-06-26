import MDAnalysis
import subprocess
import sys
import platform
def addAtoms(pdbName, pH, ripmdPath, selection, outFolder):
	if sys.platform=="linux2":	
		if platform.machine()!="x86_64":
			print "We can not found missing atoms because we can not found a PDB2PQR version for your PC architecture... omitting operation"
			return
		command=ripmdPath+"dependencies/linux/pdb2pqr-linux-bin64-2.1.1/pdb2pqr"

		#if we only work with the PDB we will perform an optimization, if not we only assign charges and vdw radii
		command+=" --ff=charmm --nodebump --noopt --chain --with-ph "+str(pH)
		pqr_name=pdbName[:-4].split("/")[-1]+".pqr" #generating pqr name
		command+=" "+pdbName+" "+outFolder+"/.temp_RIP-MD/"+pqr_name+" >> "+outFolder+"/.temp_RIP-MD/"+pqr_name+".temp" #adding this to the command, the new .temp is for not print lines in the screen and, for this way, saving compute time
		proc = subprocess.Popen([command], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True) # to call pdb2pqr in a subprocess to know if there exist any error
		out, err = proc.communicate()
		if proc.returncode==0:
				u=MDAnalysis.Universe(pdbName)
				pdb = None
				try:
					pdb=u.select_atoms(selection)
				except:
					try:
						pdb=u.selectAtoms(selection)
					except:
						
						return False
						
				up=MDAnalysis.Universe(outFolder+"/.temp_RIP-MD/"+pqr_name)
				pqr = None
				try:
					pqr=up.select_atoms(selection)
				except:
					pqr=up.selectAtoms(selection)
				#if we add some atoms, we will save the new pdb, else we will do nothing
				if len(pdb)!=len(pqr):
					print "Saving new PDB with added atoms..."
					pqr.write(outFolder+"/.temp_RIP-MD/"+(pdbName[:-4].split("/")[-1])+"_missingAtoms.pdb")
					del u
					del pdb
					del up
					del pqr
					return "good"
				else:
					print "No new atoms added... omitting operation"
		del pqr_name
		del command
		del proc
		del out
		del err
		
	else:
		print "We can not found missing atoms due a missing PDB2PQR version for your OS... omitting operation"
	
	return
	


