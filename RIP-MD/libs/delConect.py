import multiprocessing as mp
import os

def parallelDel(pdb):
	os.system("sed -i '/CONECT/d' "+pdb[1].replace(" ","\ "))



def Delete(pdbDict, nproc):
	pool=mp.Pool(processes=int(nproc))
	pool.map(parallelDel,(pdbDict.items()))

