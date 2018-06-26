import os
import glob
########################################################################
##
## library to delete a set of file with some formats
##
########################################################################

def delete(folder, Format):
	toDelete=glob.glob(folder+"/*."+Format)
	for File in toDelete:
		os.remove(File)
	
