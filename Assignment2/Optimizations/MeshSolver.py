from Matrix import Matrix
import MeshStructure
import math
from Choleski import Choleski

S_g = 0
S_ff = 0
S_fp = 0 #actually holds -1* S_fp
U_f = 0
U_p = 0
U_con_order = []

def getArea(triangle):
	#TODO i'm getting different areas in one mesh?
	n1 = triangle.node1
	n2 = triangle.node2
	n3 = triangle.node3
	
	return math.fabs( (n1.x * (n2.y - n3.y) + n2.x * (n3.y - n1.y) + n3.x * (n1.y - n2.y) )/ 2 )

# returns the inverse of a 3x3 matrix
def inverse3x3(matrix):
	a = matrix.get(1,1)
	b = matrix.get(1,2)
	c = matrix.get(1,3)
	d = matrix.get(2,1)
	e = matrix.get(2,2)
	f = matrix.get(2,3)
	g = matrix.get(3,1)
	h = matrix.get(3,2)
	k = matrix.get(3,3)
	
	determinant = a*(e*k-f*h)+b*(f*g-k*d)+c*(d*h-e*g)
	
	matrix.set(1,1,e*k-f*h)
	matrix.set(1,2,f*g-d*k)
	matrix.set(1,3,d*h-e*g)
	matrix.set(2,1,c*h-b*k)
	matrix.set(2,2,a*k-c*g)
	matrix.set(2,3,g*b-a*h)
	matrix.set(3,1,b*f-c*e)
	matrix.set(3,2,c*d-a*f)
	matrix.set(3,3,a*e-b*d)	
	
	inverse = matrix.transpose()
	
	for i in range(1,3+1):
		for j in range(1,3+1):
			inverse.set(i,j,inverse.get(i,j)*(1.0/determinant))
			
	return inverse
	
def buildLocal_S(triangle):
	area = getArea(triangle)
	
	n1 = triangle.node1
	n2 = triangle.node2
	n3 = triangle.node3
	
	# build X matrix
	X = Matrix([	[1, n1.x, n1.y],
					[1, n2.x, n2.y],
					[1, n3.x, n3.y]
				])

	# compute inverse(X)
	X_inverse = inverse3x3(X)
	
	# compute transpose(X_inverse)
	X_i_t = X_inverse.transpose()
	
	# get delta_alpha (i and j components)
	alphaI = Matrix([  	[X_i_t.get(1,2)],
						[X_i_t.get(2,2)],
						[X_i_t.get(3,2)]
					])				
	alphaJ = Matrix([	[X_i_t.get(1,3)],
						[X_i_t.get(2,3)],
						[X_i_t.get(3,3)]
					])
	
	S = alphaI.multiply(alphaI.transpose()).add(alphaJ.multiply(alphaJ.transpose()))
	
	for i in range(1,3+1):
		for j in range(1,3+1):
			S.set(i,j,S.get(i,j)*(area))
			
	return S

def findJointNodes(triangleL):
	joint=[]
	
	for k in range(0,len(triangleL)):
		k_nodes=[]
		k_nodes.append(triangleL[k].node1)
		k_nodes.append(triangleL[k].node2)
		k_nodes.append(triangleL[k].node3)
		#compare starting from k+1 as previous have already compared
		for j in range(k+1,len(triangleL)):
			j_nodes=[]
			j_nodes.append(triangleL[j].node1)
			j_nodes.append(triangleL[j].node2)
			j_nodes.append(triangleL[j].node3)
			
			#compare the two triangles
			for a in range(0,len(k_nodes)):
				for b in range(0,len(j_nodes)):
					if (k_nodes[a].x==j_nodes[b].x and k_nodes[a].y==j_nodes[b].y and k_nodes[a]!=j_nodes[b]):
						joint.append([k_nodes[a].number,j_nodes[b].number])
						
	return joint	
			
	
def buildGlobal_S(nodeH, triangleL, pList):
	global S_g
	global U_f
	global U_p
	global U_con_order
	#populate global S with zeroes
	S_g = Matrix(i=(len(triangleL)*3),j=(len(triangleL)*3))
	for i in range(1,S_g.rows+1):
		for j in range(1,S_g.columns+1):
			S_g.set(i,j,0)
	
	#build local S for each triangle and place it in global S
	i = 1 #index where we'll be placing the local S into global S
	j = 1
	for k in range(0,len(triangleL)):
		S_l = buildLocal_S(triangleL[k])
		#add it to S_g in appropriate place
		for a in range(1,S_l.rows+1):
			for b in range(1, S_l.columns+1):
				S_g.set(i,j,S_l.get(a,b))
				j+=1
			i+=1
			j-=3
		j+=3
	
	#No need for findJoint in the case of our meshes as nodes are already joined	
	#joint = findJointNodes(triangleL)
	
	#build a hash from id(node) -> node
	nodeID = {}
	for i in range(1,len(nodeH)+1):
		nodeID[id(nodeH[i])] = nodeH[i]
	
	#build U disjoint with node ids (rather than voltages)
	U_dis = Matrix(i=(len(triangleL)*3),j=1)
	row = 1
	for k in range(0,len(triangleL)):
		U_dis.set(row,1,id(triangleL[k].node1))
		U_dis.set(row+1,1,id(triangleL[k].node2))
		U_dis.set(row+2,1,id(triangleL[k].node3))
		row+=3
		
	
	#build U conjoint with node ids (and build C along the way) 
	C = Matrix(i=U_dis.rows,j=len(nodeH))
	row_C = 1
	U_con = Matrix(i=len(nodeH),j=1)
	row_U_c = 1 #row in conjoint where we place our next element
	for row_d in range(1,U_dis.rows+1):
		node_id = U_dis.get(row_d,1)
		add = True
		#check if it's already in U_con
		for j in range(row_U_c-1,0,-1):
			if (U_con.get(j,1) == node_id):
				add = False
				break
		if (add):
			#add it to U_con since it hasn't been added yet
			U_con.set(row_U_c,1,node_id)
			C.set(row_C,row_U_c,1.0)
			row_U_c+=1			
			row_C+=1
		else:
			C.set(row_C,j,1.0)
			row_C+=1
	
	#rearrange	U_con and C (to put free up, and predefined down)
	#method used is 2 pointers moving in opposite directions on U_con
	#the pointer moving up looks for a free node, and the pointer moving down looks for predef node
	#swap both, stop process when cross
	down = 1 #index moving down
	up = U_con.rows #index moving up
	while (down < up):
		while (down < up and nodeID[U_con.get(down,1)].number not in pList):
			down+=1
		if (nodeID[U_con.get(down,1)].number not in pList):
			break #nothing more to switch around
			
		while (up > down and nodeID[U_con.get(up,1)].number in pList): #looking for free from bottom
			up-=1
		
		if (nodeID[U_con.get(up,1)].number in pList):
			break #nothing more to switch around
		
		#do the swap for U_con
		tmp = U_con.get(up,1)
		U_con.set(up,1,U_con.get(down,1))
		U_con.set(down,1,tmp)
		
		#do the swap for C (swap columns)
		for i in range(1,C.rows+1):
			tmp = C.get(i,up)
			C.set(i,up,C.get(i,down))
			C.set(i,down,tmp)
	
	#convert U_con to actual voltage values (but keep track of ordering)
	for i in range(1,U_con.rows+1):
		U_con_order.append(nodeID[U_con.get(i,1)].number)
		U_con.set(i,1,nodeID[U_con.get(i,1)].voltage)
	
	#build U_p and U_f
	U_f = Matrix(i=down-1,j=1)
	U_p = Matrix(i=U_con.rows-U_f.rows,j=1)
	for i in range(1,U_f.rows+1):
		U_f.set(i,1,U_con.get(i,1))
	for i in range(1,U_p.rows+1):
		U_p.set(i,1,U_con.get(i+down-1,1))
	
	#finally compute S_g
	S_g = (C.transpose().multiply(S_g)).multiply(C)

def reduced():
	global S_ff
	global S_fp
	
	S_ff = Matrix(i=(U_f.rows),j=(U_f.rows))
	S_fp = Matrix(i=(U_f.rows),j=(U_p.rows))
	for i in range(1,S_ff.rows+1):
		for j in range(1,S_g.columns+1):
			if (j <= U_f.rows):
				S_ff.set(i,j,S_g.get(i,j))
			else:
				S_fp.set(i,j-U_f.rows,-1*S_g.get(i,j))

def printResults(nodeH,pList,filename):
	f = open(filename,"w")
	
	#put results in appropriate nodes in nodeH
	for i in range(1,U_f.rows+1):
		nodeH[U_con_order[i-1]].voltage = U_f.get(i,1)
		
	f.write("node,x,y,potential,predefined\n")
	for i in range(1,len(nodeH)+1):
		f.write(str(nodeH[i].number)+","+str(nodeH[i].x)+","+str(nodeH[i].y)+","+str(nodeH[i].voltage))
		if ( nodeH[i].number in pList):
			f.write(",y\n")
		else:
			f.write("\n")
		
	f.close()
				
#t1=MeshStructure.Triangle(MeshStructure.Node(1,0,0.02,0.0),MeshStructure.Node(2,0,0,0.0),MeshStructure.Node(3,0.02,0,0.0))
#t2=MeshStructure.Triangle(MeshStructure.Node(4,0.02,0.02,0.0),MeshStructure.Node(5,0,0.02,0.0),MeshStructure.Node(6,0.02,0,0.0))
#buildGlobal_S([t1,t2])
#MeshStructure.buildStructs('squareCoax3_D.msh')
import sys
if (len(sys.argv) > 1):
	filename = sys.argv[1]
	MeshStructure.buildStructs(filename+".msh")
	buildGlobal_S(MeshStructure.nodeH, MeshStructure.triangleL, MeshStructure.pList)
	reduced()
	U_f = Choleski(S_ff,S_fp.multiply(U_p))
	printResults(MeshStructure.nodeH,MeshStructure.pList,filename+".csv")
else:
	print "please specify a filename without extension"





