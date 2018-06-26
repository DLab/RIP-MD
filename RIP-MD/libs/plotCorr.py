import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis
import sys

####################################################
##
## Main function
##
####################################################
def Pearson(pdb, output, file,fromRIPMD=False):

	matrix = np.loadtxt(file,skiprows=0)

	##################################################
	##
	## looking for chain and IDs to generate
	## plot axis
	##
	##################################################
	chains=[[],[]]
	u=MDAnalysis.Universe(pdb)
	sel=u.select_atoms("protein").segments
	total=0
	for segment in sel:
		chains[0].append(segment.segid)
		chains[1].append(len(segment.atoms.residues))
		total+=len(segment.atoms.residues)
	x_axis_lenght=((-1*(total))/8)
	fig, ax = plt.subplots()	
	#adding axis
	x1=x_axis_lenght
	x2=x_axis_lenght
	x3=((-1*(total))/5)
	y1=0
	y2=0
	x3=((-1*(total))/5)
	y3=0
	X1=0
	Y1=total+(total*0.06)
	X2=0
	Y2=total+(total*0.06)
	X3=0
	Y3=total+(total*0.1)
	text=""
	for i in range (len(chains[0])):
		#adding y axis
		y1=0
		if i>0:
			for j in range(i):
				y1+=chains[1][j]
		
		y2=y1+chains[1][i]
		
		y3=(y1+y2)/2.
		text=chains[0][i]
		ax.annotate('', xy=(x1, y1),xytext=(x2,y2), 	#draws an arrow from one set of coordinates to the other
		arrowprops=dict(arrowstyle='<->'),				#sets style of arrow
		annotation_clip=False)							#This enables the arrow to be outside of the plot	
		
		ax.annotate(text,xy=(0,0),xytext=(x3,y3),               #Adds another annotation for the text
		annotation_clip=False)

		#adding x axis

		X1=0
		if i>0:
			for j in range(i):
				X1+=chains[1][j]
		X2=X1+chains[1][i]
		X3=(X1+X2)/2.
		ax.annotate('', xy=(X1, Y1),xytext=(X2,Y2), 	#draws an arrow from one set of coordinates to the other
		arrowprops=dict(arrowstyle='<->'),				#sets style of arrow
		annotation_clip=False)							#This enables the arrow to be outside of the plot
		
		ax.annotate(text,xy=(0,0),xytext=(X3,Y3),               #Adds another annotation for the text
		annotation_clip=False)				
	plt.imshow(matrix, interpolation="none")
	plt.clim(-1,1)

	plt.colorbar()


	if fromRIPMD==False:
		#matplotlib.use('Agg')
		plt.show()		
	else:
		aux=file.split("/")
		matplotlib.use('Agg')
		plt.savefig(output+"/RIP-MD_Results/Correlations/"+aux[-1]+".png")

if __name__ == "__main__":
	Dict={"All":"total","C Alpha":"ca", "H Bonds":"hb", "Salt Bridges":"salt","Cation - pi interaction":"cation-pi", "pi - pi interaction":"pi-pi","Arg - Arg interaction":"arg-arg", "Coulomb interaction":"coulomb", "VdW contacts":"vdw"}
	file1=sys.argv[1]
	file2=sys.argv[2]
	folder=sys.argv[3]
	interaction1=Dict[file1]
	interaction2=Dict[file2]
	
	Pearson(folder+"/RIP-MD_Results/representative.pdb", None, folder+"/RIP-MD_Results/Correlations/"+interaction1+"_"+interaction2)
#Pearson("/home/scontreras/Escritorio/result_MD2/RIP-MD_Results/representative.pdb", None, "/home/scontreras/Escritorio/result_MD2/RIP-MD_Results/Correlations/ca_ca")		
#Pearson("/home/scontreras/Escritorio/result_gjc/RIP-MD_Results/representative.pdb", None, "/home/scontreras/Escritorio/result_gjc/RIP-MD_Results/Correlations/ca_ca")
