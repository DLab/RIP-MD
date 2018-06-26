import argparse   # arguments parser
import time
import os
import sys
import shutil
import MDAnalysis as mda
from platform import platform

## appending path to my python path ##
sys.path.append(str(os.path.realpath(__file__))[:-8]+"/libs") # path where this script is located
from selection import save_selection
import missingAtoms
import delConect
import checkPDB
import extractFrames as ef
import delete_hetatm as dh  #..to clean pdbs
import pdb2graph as p2g
import nodes
import Calpha as ca
import hbond as hb
import saltBridge as sb
import disulphide as ss
import cationPi as cp
import pi_pi
import arg_arg
import coulomb
import vdw
import graph

import correlation
import plotCorr
ripmdPath=os.path.realpath(__file__)[:-7] #path to RIP-MD
if len(sys.argv)<=1:
	print "use --help  option for usage information... exiting"
	exit()

############################################
##
## definition of arguments
##
############################################

#input
parser = argparse.ArgumentParser()

parser.add_argument("-V", "--version", help="Print RIP-MD version and exit",action="store_true")

parser.add_argument("-D", "--dcd", help="Route to the DCD file")
parser.add_argument("-S", "--psf", help="Route to the PSF file")
parser.add_argument("-P", "--pdb", help="Route to PDB file")
parser.add_argument("-sel","--selection",help="Protein selection in CHARMM like format (as defined in MDAnalysis) to perform RINs. Default: \"protein\"", default="protein")

#output route
parser.add_argument("-O","--output", help="Output path to save results")

#number of processors to use
parser.add_argument("-N","--nproc",help="Number of processors to use. Default: 1", default="1")

#general options for PDB2PQR and others
parser.add_argument("-ma","--missing_atoms", help="Only for single PDB: Add missing atoms (including H atoms) though PDB2PQR software", action="store_true")
parser.add_argument("-pH","--pH", help="Only for single PDB: pH for the compute missing atoms step. Default: 7.0", default="7.0")
parser.add_argument("-ff","--force_field", help="path to CHARMM forcefield to use. By default we will use our local copy of the CHARMM forcefield ", default=ripmdPath+"dat/top_all22_prot.rtf")
parser.add_argument("-pf","--parameter_file", help="path to CHARMM parameter file. By default we will use our local copy of the CHARMM parameter file", default=ripmdPath+"dat/par_all22_prot.prm")

#options for Calpha
parser.add_argument("-ca","--calpha",help="Calculate distances between 2 alpha carbons. options are [-Cad | --Ca_dist]", action="store_true")
parser.add_argument("-cad","--ca_dist",help="Distance cutoff in angstromto use when creating alpha carbons network (default  8)", default="8.0")

#options for HBond
parser.add_argument("-H","--hbond",help="Calculate Hbond network of protein. options for H Bond are [-hd | --h_dist] and [-ha | --h_angle]", action="store_true")
parser.add_argument("-hd","--h_dist",help="Distance cutoff in angstroms to use when creating hbond network (default 3)", default="3")
parser.add_argument("-ha","--h_angle",help="Angle cutoff in degrees to use when creating hbond network (default 120)", default="120")

#option for Salt Bridges
parser.add_argument("-s","--salt",help="calculate Salt bridge network between residues. options for Salt bridge network are [-sd | --s_dist]",  action="store_true")
parser.add_argument("-sd","--s_distance",help="Distance cutoff in angstromto use when creating salt bridge network (default 6)", default="6.0")

#option for disulfide bond
parser.add_argument("-ss","--disulfide",help="Calculate disulfide bond between 2 CYS. options for disulfide network are [-ssd | --ss_distance] and [-ssa | --ss_angle]", action="store_true")
parser.add_argument("-ssd","--ss_distance",help="Distance cutoff in angstrom to use when creating disulfide bridge network (default 3)", default="3.0")
parser.add_argument("-ssa","--ss_angle",help="Angles ranges between vectors formed by C-S atoms in degrees to consider disulfide bond. Accepted format is [X,Y] where X is the first angle and Y is the second angle (default [60,90])", default="[60,90]")

#option for cation-pi
parser.add_argument("-cp","--cation_pi",help="Calculate cation-pi interaction. options are [-cpd | --cation_pi_distance], [-cpa1 | --cp_angle1], [-cpa2 | --cp_angle2], [-cpa3 | --cp_angle3] and [-cpa4 | --cp_angle4]", action="store_true")
parser.add_argument("-cpd","--cation_pi_distance",help="Distance cutoff in angstrom to use when creating cation-pi network (default 7)", default="7.0")
parser.add_argument("-cpa1","--cp_angle1",help="Angles in degrees ranges between the cation mass center and the vector normal to the aromatic ring plane in format [X1,Y1,X2,Y2] (X1<alpha<Y1 or X2<alpha<Y2) whithout spaces, where X is the first angle and Y is the second angle. Default [0,60,120,180]", default="[0,60,120,180]")
parser.add_argument("-cpa2","--cp_angle2",help="Angles in degrees ranges that define the planar spatial configuration in format [X1,Y1,X2,Y2] (X1<alpha<Y1 or X2<alpha<Y2) whithout spaces, where X is the first angle and Y is the second angle. Default [0,30,150,180]", default="[0,30,150,180]")
parser.add_argument("-cpa3","--cp_angle3",help="Angles in degrees ranges that define the oblique spatial configuration in format [X1,Y1,X2,Y2] (X1<alpha<Y1 or X2<alpha<Y2) whithout spaces, where X is the first angle and Y is the second angle. Default [30,60,120,150", default="[30,60,120,150]")
parser.add_argument("-cpa4","--cp_angle4",help="Angles in degrees ranges that define the orthogonal spatial configuration in format [X1,Y1] (X1<alpha<Y1) whithout spaces, where X is the first angle and Y is the second angle. Default [60,120]", default="[60,120]")

#option for pi-pi
parser.add_argument("-pp","--pi_pi",help="Calculate pi-pi interaction. options are [-ppd | --pp_distance]", action="store_true")
parser.add_argument("-ppd","--pi_pi_distance",help="Distance cutoff in angstroms to use when creating pi-pi network (default 6)", default="6.0")

#option for arginine-arginine
parser.add_argument("-rr","--arg_arg",help="Calculate Arginine-Arginine interaction. options are [-rrd | --arg_arg_distance]", action="store_true")
parser.add_argument("-rrd","--arg_arg_distance",help="Distance cutoff in angstroms to use when creating arginine-arginine network (default 5)", default="5.0")

#option for vdw
parser.add_argument("-v","--vdw",help="Calculate van der Waals interactions. options are [-vd | --vdw_distance]", action="store_true")
parser.add_argument("-ve","--vdw_excluded",help="Covalent Bond to exclude vdw calculus (Distance cutoff to create vdW contact network. Default 3 (that represent 3 bonds or 1-4 potential)", default="3")
parser.add_argument("-vr","--vdw_range",help="maximum and minimum distance in angstroms between vdw radii of two atoms in format [X,Y] (default [-0.1,3] where -0.1 are two atoms overlapping)", default="[-0.1,3]")

#option for coulomb potential
parser.add_argument("-c","--coulomb",help="Calculate Coulomb potential between charged groups of each aminoacids defined on the CHARMM force field", action="store_true")
parser.add_argument("-RF","--reaction_field",help="Add a reaction field formulation to evaluate the long distance  electrostatic interaction", action="store_true")
parser.add_argument("-Ecs","--simulated_permittivity",help="Value for the simulated permittivity. Default value: 1 ", default="1")
parser.add_argument("-Erf","--RF_permittivity",help="Value for the reaction field permittivity. Default value: 80 ", default="82")
parser.add_argument("-T","--temperature", help="Temperature used to compute thermal noise cutoff. Default 293.5", default="298.5")
parser.add_argument("-Krf","--inverse_debye_screening", help="Inverse Bebye screening length. Default 0", default="0")
parser.add_argument("-Rrf","--distance_threshold",help="Coulomb distance threshold to define interaction. Default value 12", default="12")
parser.add_argument("-ce","--coulomb_excluded",help="Covalent Bond to exclude coulomb calculus (Distance cutoff to create coulombic network. Default 3 (that represent 3 bonds or 1-4 potential)", default="3")


#options for consensus graph
parser.add_argument("-t","--time", help="Percentage of time where interactions must be present (default 75)", default="75")
parser.add_argument("-gf","--gformat", help="Format for the output graphs. Options are:(1) edgelist || (2) GEXF || (3) GML || (4) GraphML || (5) Pajek. Default option: GML", default="GML")

#options for frame windows to get correlation between aminoacid interactions
parser.add_argument("-pc","--pearson_corr",help="Calculate Pearson Correlation between residues interactions", action="store_true")
parser.add_argument("-p","--plot_pearson", help="generate plots for computed correlations", action="store_true")
#options to save frame
parser.add_argument("-sf","--separation_frame", help="frame separation to compute interactions. Default:0",default="1")
parser.add_argument("-fs","--frame_start", help="frame number to start frame extraction. Default:0",default="0")
parser.add_argument("-fe","--frame_end", help="frame number to finich frame extraction. Default: last frame (value -1)",default="-1")
parser.add_argument("-rf","--reference_frame",help="Reference frame to display interactions. If your reference frame does not exist, the last frame will be used. Default Frame: 0", default="0")


args = parser.parse_args()

if args.version:
	print "RIP-MD Version 2.0.0"
	exit()

#######################################################
##
## error options
##
#######################################################

if not os.path.exists(args.output): #if path does not exist
	print "Output folder does not exist, creating..."
	try:
		os.mkdir(args.output)
	except:
		print "Can not create directory "+args.output+" ... exiting"
		exit()
sys.stdout = open(args.output+"/output.log", 'w')
sys.stderr = open(args.output+"/output.log", 'a')


if args.dcd and args.psf and args.pdb:
	print "Error: you can not provide a dcd, psf and a pdb file as input"
	exit()

elif args.dcd and args.pdb:
	print "Error: you can not provide a dcd and a pdb file as input"
	exit()

elif args.dcd and args.pdb:
	print "Error: you can not provide a dcd and a pdb file as input"
	exit()

elif args.dcd and not args.psf and not args.pdb:
	print "Error: you can not provide only a dcd file as input"
	exit()

elif not args.dcd and args.psf and not args.pdb:
	print "Error: you can not provide only a psf file input"
	exit()

else:
	startMssg= "Starting RIP-MD using command python"
	for Args in sys.argv:
		startMssg+=" "+str(Args)
	print startMssg
	start_time = time.time()
	print "Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+""
	print "Creating folders for temporal and results files..."
	try:
		os.mkdir(args.output+"/.temp_RIP-MD")
		print "\tTemporal folder created..."
	except:
		print "\tTermporal folder already exist... omitting..."
	try:
		if args.output:
			os.mkdir(args.output+"/RIP-MD_Results")
		print "\tResults folder created..."
	except:
		try:
			shutil.rmtree(args.output+"/RIP-MD_Results/Attrs")
		except:
			pass
		try:
			shutil.rmtree(args.output+"/RIP-MD_Results/Edges")
		except:
			pass
		try:
			shutil.rmtree(args.output+"/RIP-MD_Results/Graphs")
		except:
			pass
		#os.mkdir(args.output+"/RIP-MD_Results")
		print "\tResults folder created..."
	try:
		os.mkdir(args.output+"/RIP-MD_Results/Edges")
		print "\tEdges folder created..."
	except:
		print "\tEdges folder already exist... omitting..."
	try:
		os.mkdir(args.output+"/RIP-MD_Results/Graphs")
		print "\tGraphs folder created..."
	except:
		print "\tGraphs folder already exist... omitting..."
	#######################################################
	##
	## global variables
	##
	#######################################################
	pdbNames={} # a dictionary of frame:pdbName
	if args.pdb:

		"""/////////////////////////////////////////////////////////////////
		""
		"" adding missing atoms
		""
		/////////////////////////////////////////////////////////////////"""
		if args.pdb:
			pdbNames["frame_0"] = args.pdb
			if args.missing_atoms and not args.psf:
				print "Checking and adding missing atoms with PDB2PQR software..."
				MA = missingAtoms.addAtoms(pdbNames.items()[0][1],args.pH,ripmdPath, args.selection, args.output)
				if MA =="good":
					pdbNames["frame_0"] = args.output+"/.temp_RIP-MD/"+(pdbNames["frame_0"][:-4].split("/")[-1])+"_missingAtoms.pdb"
				if MA == False:
					print "Error adding atoms for the current selection, exiting (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."
					exit()
				print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"
				
			##########################################
			##
			## Now we will create a PDB with selection
			##
			##########################################
				
			print "Creating a copy of your PDB with current selection: "+args.selection
			save = save_selection(args.selection, pdbNames["frame_0"], args.output+"/.temp_RIP-MD/"+str(args.pdb[0:-4].split("/")[-1])+"_selection.pdb")
			if save != False:
				pdbNames["frame_0"] = save
				print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"
			else:
				print "Can not create selection, exiting (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."
				exit()			

	elif args.dcd and args.psf:
		try:
			os.mkdir(args.output+"/RIP-MD_Results/Correlations")
			print "\tCorrelations folder created..."
		except:
			print "\tCorrelations folder already exist... omitting..."
		########################################################################
		##
		## For MD, we will extract frames and then we will complete the "pdb
		## dictionary" with all names
		##
		########################################################################
		print "Extracting Frames according to selection (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."
		pdbNames=ef.DCD(args.psf, args.dcd, args.output, args.separation_frame, args.frame_start, args.frame_end, args.nproc, args.selection)		
		print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"
		
	
	"""/////////////////////////////////////////////////////////////////
	""
	"" deleting conect lines in pdbs
	""
	/////////////////////////////////////////////////////////////////"""

	#we will delete these lines due mdanalysis write conects with numbers very close in the text
	# so it cause errors when it it read pdb
	print "Deleting CONECT lines from PDBs..."
	delConect.Delete(pdbNames, args.nproc)
	print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"

	"""/////////////////////////////////////////////////////////////////
	""
	"" Checking for repeated residues
	""
	/////////////////////////////////////////////////////////////////"""
	if args.psf==None:
		print "Checking for problems with your structure..."
		check=checkPDB.check(args.output, pdbNames, args.nproc)
		if check==1: ##if exist repeated residues we exit from RIPMD
			print " (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"
			exit()
		print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"	

	"""/////////////////////////////////////////////////////////////////
	""
	"" at this point we have all necessary files to extract information,
	"" so we will take the first pdb of the list and we will save amino
	"" acids. then, for each pdb we will save a file with attributes
	"" and edges
	""
	/////////////////////////////////////////////////////////////////"""

	print "Parsing Nodes (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."
	Nodes=nodes.write_nodes(args.output,pdbNames["frame_0"])# Nodes var will have id of node and the node (only take the first frame)
	print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"

	## Now we are going to compute DSSP
	print "Writting attributes (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."
	#outDSSP will give us information about if DSSP finished with any error
	outDSSP=nodes.secStructure(pdbNames, Nodes, args.output, args.nproc, ripmdPath)
	print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"
	
	
	"""///////////////////////////////////////
	""
	"" Solving structure topology
	""
	///////////////////////////////////////"""
	topology=None
	if args.vdw or args.coulomb:
		print "Trying to solve structure topology. This could take several time (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."
		topology=p2g.convert(pdbNames["frame_0"], args.psf, args.force_field,args.parameter_file, args.output)
		print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"

	"""///////////////////////////////////////
	""
	"" calculating distance between Calphas
	""
	//////////////////////////////////////"""
	if args.calpha:
		print "Computing contact map based on C alpha (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."
		ca.parse(pdbNames, args.nproc, args.output, Nodes, float(args.ca_dist))
		print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"

	"""///////////////////////////////////////
	""
	"" calculating HBond
	""
	//////////////////////////////////////"""
	if args.hbond: ##if hbond option is true
		print "Computing HBonds (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."
		hb.parse(pdbNames,args.output, args.nproc, Nodes, ripmdPath, args.h_dist, args.h_angle) #we will call HB function and we will parse the universe, distance and angle
		print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"

	"""///////////////////////////////////////
	""
	"" Calculating Salt Bridges
	""
	//////////////////////////////////////"""
	if args.salt:
		print "Computing Salt bridges (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."
		sb.parse(pdbNames, args.output, args.nproc, Nodes, float(args.s_distance)) 	
		print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"

	"""//////////////////////////////////////
	""
	"" calculating disulphide bond
	""
	//////////////////////////////////////"""
	if args.disulfide:
		if args.psf:
			print "RIP-MD will not determine disulfide bonds for your input because they will be kept through all the simulation due to definition in the PSF file... omitting. (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"
		else:
			print "Computing disulfide bonds (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."	
			ss.parse(pdbNames,args.output, args.nproc, Nodes, float(args.ss_distance), args.ss_angle)
			print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"

	"""//////////////////////////////////////
	""
	"" calculating cation pi interactions
	""
	//////////////////////////////////////"""
	if args.cation_pi:
		print "Computing Cation-pi interaction (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."
		cp.parse(pdbNames,args.psf, args.output, args.nproc, Nodes, float(args.cation_pi_distance), args.cp_angle1, args.cp_angle2, args.cp_angle3, args.cp_angle4)
		print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"

	"""//////////////////////////////////////
	""
	"" calculating pi-pi interactions
	""
	//////////////////////////////////////"""
	if args.pi_pi:
		print "Computing pi-pi interaction (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."		
		pi_pi.parse(pdbNames, args.output, args.nproc, Nodes, float(args.pi_pi_distance))
		print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"
	"""//////////////////////////////////////
	""
	"" calculating arg-arg interactions
	""
	//////////////////////////////////////"""	
	if args.arg_arg:
		print "Computing ARG-ARG interaction (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."		
		arg_arg.parse(pdbNames, args.output, args.nproc, Nodes, float(args.arg_arg_distance))
		print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"		

	"""//////////////////////////////////////
	""
	"" calculating vdW contacts
	""
	//////////////////////////////////////"""	
	if args.vdw:
		print "Computing VdW contacts (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."	
		vdw.parse(args, pdbNames, args.output, args.nproc, Nodes, args.vdw_excluded, args.vdw_range, topology)
		print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"

	"""//////////////////////////////////////
	""
	"" calculating coulomb interactions
	""
	//////////////////////////////////////"""	
	if args.coulomb:
		print "Computing Coulomb interactions (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."	
		coulomb.parse(pdbNames, args.output, args.nproc, Nodes, topology, float(args.simulated_permittivity), float(args.RF_permittivity), float(args.temperature), float(args.inverse_debye_screening), float(args.distance_threshold), float(args.coulomb_excluded), args.reaction_field, args.selection)
		print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"

	"""/////////////////////////////////////
	""
	"" doing graphs
	""
	////////////////////////////////////"""
	print "Generating graphs (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."	
	graph.doGraph(pdbNames, args.output, args.gformat, args, args.nproc, args.time) # args is to know what kind of interactions were calculated
	print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"

	"""/////////////////////////////////////
	""
	"" doing pearson correlation
	""
	////////////////////////////////////"""
	if len(pdbNames)>1:
		if args.pearson_corr:
			print "Computing correlations (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."
			corr=correlation.Pearson(pdbNames, args.output, Nodes, args.nproc, args) # args is to know what kind of interactions were computed
			print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"
			if args.plot_pearson:
				print "Plotting correlations (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."
				for file in corr:
					plotCorr.Pearson(pdbNames["frame_0"], args.output, file,True) # args is to know what kind of interactions were computed
				print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"

	#now we will look for the last frame and we move to result folder, then we will delete .temp folder
	lastPDB=None
	if args.pdb:
		lastPDB=str(pdbNames["frame_"+str(len(pdbNames)-1)])
		shutil.move(lastPDB,args.output+"/RIP-MD_Results/representative.pdb")
	else:
		print "Copying the representative frame (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")..."

		try:
			lastPDB=str(pdbNames["frame_"+str(args.reference_frame)])
		except:
			print "Not possible to use your reference frame, using last frame...\n"
			lastPDB=str(pdbNames["frame_"+str(len(pdbNames)-1)])
		shutil.move(lastPDB,args.output+"/RIP-MD_Results/representative.pdb")
		print "Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")"
	
	try:
		shutil.rmtree(args.output+"/.temp_RIP-MD")
		print "Deleting temporal files"
	except:
		pass

	file=open(args.output+"/README","w")
	file.write("Result Description:\n\n")
	file.write("In the output folder there is a folder called \"RIP-MD_Results\" and a log file (output.log) where it is possible to see the status of RIP-MD.\n\nInside the RIP-MD_Results folther there are several other folders.\n\n")
	file.write("Inside the Attrs folder, there are files describing the attributes for the loaded PDB or for each frame of the MD.")
	file.write("Inside the Edges folder, there are files describing the interactions formed for the loaded PDB or for each frame of the MD.")
	file.write("Inside the Graphs folder, there are files with representation of each type of interaction between 2 or more aminoacids named after the interaction type.")
	file.write("\n\n\nFile format description:\n\n")
	file.write(args.gformat+" files: Network files for easy loading and visualization in Cytoscape 3.X. These files contain each aminoacid represented as a node and the interactions (and disulphide bonds) represented as edges.\n")
	file.write("attr files: A simple text file in tsv format with node attributes for the respective frame.\n")
	file.write("edges files: A simple text file in tsv format where edge and interaction properties are described for each frame.\n")
	file.write("pearson files: A simple text file in tsv format with pearson's correlation between nodes.\n")
	file.write("png file: A image to visualize pearson's correlation between nodes.\n")
	file.close()

print "END"
