import math
import numpy as np
def angle(v1,v2,v3): #v1 and v2 are two points in the space, v3 is a normal vector
	#first we determine a auxiliar vector (v2-v1)
	v_aux=v2-v1
	
	#now we determine the escalar vector
	scalar_prod=(v_aux[0]*v3[0])+(v_aux[1]*v3[1])+(v_aux[2]*v3[2])
	#now we determine the magnitude
	m_vaux= math.sqrt((v_aux[0]**2)+(v_aux[1]**2)+(v_aux[2]**2))
	m_v3=math.sqrt((v3[0]**2)+(v3[1]**2)+(v3[2]**2))
	
	#and now we will find the angle
	angle=1/(math.cos(scalar_prod/(m_vaux*m_v3)))
	return angle*57.2957795
	
def angle_polar_groups(v1,v2,v3,v4): #each v is a point in the space
	#first we determine a auxiliar vector (v2-v1)

	v_aux=v2-v1
	v_aux2=v4-v3
	
	#now we determine the escalar vector
	scalar_prod=(v_aux[0]*v_aux2[0])+(v_aux[1]*v_aux2[1])+(v_aux[2]*v_aux2[2])
	#now we determine the magnitude
	m_vaux= math.sqrt((v_aux[0]**2)+(v_aux[1]**2)+(v_aux[2]**2))
	m_vaux2=math.sqrt((v_aux2[0]**2)+(v_aux2[1]**2)+(v_aux2[2]**2))
	
	#and now we will find the angle
	angle=1/(math.cos(scalar_prod/(m_vaux*m_vaux2)))
	return angle*57.2957795	
	
def angle_twoVectors(v1,v2):
	
	#now we determine the escalar vector
	scalar_prod=(v1[0]*v2[0])+(v1[1]*v2[1])+(v1[2]*v2[2])
	#now we determine the magnitude
	m_v1= math.sqrt((v1[0]**2)+(v1[1]**2)+(v1[2]**2))
	m_v2=math.sqrt((v2[0]**2)+(v2[1]**2)+(v2[2]**2))
	
	#and now we will find the angle
	angle=1/(math.cos(scalar_prod/(m_v1*m_v2)))
	angle=angle*57.2957795	
	if angle<0:
		angle+=180
	return angle

def dihedral(v1,v2,v3,v4):
	####################################################################
	## 
	## to computhe dihedral angel, we will compute two cross product
	## (normal vector) and then we will compute angle between these
	## normal vector
	##
	####################################################################
	
	## this function compute dihedral angle between atoms in this form
	##  a--b
	##      \
	##       c--d
	## so the first normal vector is defined by the planes (a-b)x(c-b)
	## and the second one is defined by (b-c)x(d-c)
	b0 = -1.0*(v2 - v1)
	b1 = v3 - v2
	b2 = v4 - v3
	# normalize b1 so that it does not influence magnitude of vector
	# rejections that come next
	b1 /= np.linalg.norm(b1)
	# vector rejections
	# v = projection of b0 onto plane perpendicular to b1
	#   = b0 minus component that aligns with b1
	# w = projection of b2 onto plane perpendicular to b1
	#   = b2 minus component that aligns with b1
	v = b0 - np.dot(b0, b1)*b1
	w = b2 - np.dot(b2, b1)*b1
	# angle between v and w in a plane is the torsion angle
	# v and w may not be normalized but that's fine since tan is y/x
	x = np.dot(v, w)
	y = np.dot(np.cross(b1, v), w)
	if np.degrees(np.arctan2(y, x)) < 0:
		return np.degrees(np.arctan2(y, x)) + 180
	else:
		return np.degrees(np.arctan2(y, x))
