import math
########################################################################
##
## function that take 2 arrays and calculate euclidian distance
##
########################################################################
def euclidean(array1, array2):
	x1=float(array1[0])
	y1=float(array1[1])
	z1=float(array1[2])
	
	x2=float(array2[0])
	y2=float(array2[1])
	z2=float(array2[2])
	
	distance=math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
	return distance
