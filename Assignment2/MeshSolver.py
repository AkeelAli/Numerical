from Matrix import Matrix
import MeshStructure
import math

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
	
def buildGlobal_S(triangleL):
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
		
	print S_g

t1=MeshStructure.Triangle(MeshStructure.Node(1,0,0.02,0.0),MeshStructure.Node(2,0,0,0.0),MeshStructure.Node(3,0.02,0,0.0))
t2=MeshStructure.Triangle(MeshStructure.Node(4,0.02,0.02,0.0),MeshStructure.Node(5,0,0.02,0.0),MeshStructure.Node(6,0.02,0,0.0))
buildGlobal_S([t1,t2])
#MeshStructure.buildStructs('squareCoax3_D.msh')
#buildGlobal_S(MeshStructure.triangleL)

