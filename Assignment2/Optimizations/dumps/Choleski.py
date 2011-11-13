import math
import copy
from time import clock
from Matrix import Matrix

#Solves and prints x in Ax=b given matrices A and b
def Choleski(A,b,halfBandwidth=None):
	#if no halfBandwidth is provided, then it equals columns=rows
	if (halfBandwidth is None):
		halfBandwidth=A.columns
		
	originalA=copy.deepcopy(A)
	originalb=copy.deepcopy(b)
	
	n=A.rows
	
	#checks to see that A is square, and b matches A
	assert n==A.columns
	assert n==b.rows
	
	#TODO checks to see if A is symmetric?
	
	startTime=clock()
	
	#ELIMINATION
	for j in range(1,n+1):
		if (A.get(j,j)<=0): 
			print "Choleski Error: Passed A is not real, symmetric or positive definite"
			return -1
		
		A.set(j,j,math.sqrt(A.get(j,j)))
		b.set(j,1,b.get(j,1)/A.get(j,j))

		#added this to set the upper part of A (L) to 0 [optional in the algorithm]
		for i in range(1,j):
			A.set(i,j,0)
		
		for i in range(j+1,halfBandwidth+1):
			A.set(i,j,A.get(i,j)/A.get(j,j))
			b.set(i,1,b.get(i,1)-(A.get(i,j)*b.get(j,1)))
			
			for k in range(j+1,i+1):
				A.set(i,k,A.get(i,k)-(A.get(i,j)*A.get(k,j)))
				
	#BACK SUBSTITUTION
	x=Matrix(i=n,j=1) #create empty matrix of specified dimensions
	
	#counting backwards starting from n
	for i in range(n,0,-1):
		#compute the summation
		sum=0
		for j in range(i+1,halfBandwidth+1):
			sum+=A.get(j,i)*x.get(j,1)
		
		x.set(i,1,(b.get(i,1)-sum)/A.get(i,i))
	
	elapsedTime=clock()-startTime
	
	# print "=============================================="
	# print "START Choleski Function Output						"
	# print "=============================================="
	# print ""
	# print "Given:"
	# print "-------------------------"
	# print "A ="
	# print originalA
	# print "b ="
	# print originalb
	# print ""
	# print "-------------------------"
	# print "Solution"
	# print "-------------------------"
	# print "x ="
	# print x
	# print "=============================================="
	# print "END Choleski Function Output						"
	# print "=============================================="
	print ""
	
	return x
	#return elapsedTime

#gets the half bandwidth of a symmetric matrix A
#assumes A symmetric
def getHalfBandwidth(A):
	if (A.rows!=A.columns):
		print "Matrix provided to getHalfBandwidth is not symmetric"
		return -1
		
	hb=0
	for i in range(1,A.rows+1):
		for j in range(1,i+1):
			if (A.get(i,j)!=0):
				#set hb=max(hb,i-j+1)
				if (hb < (i-j+1)):
					hb=(i-j+1)
					break;
					
	return hb

#takes lower matrix L, generates real, symmetric and positive definite A (A=LL_T)
#gets b by multiplying x by A (b=Ax)
#calls Choleski and prints expected vs obtained results
#CONDITION: x must be vector of length=rows of lower matrix L
def testCholeski(L,x):
	L_T=L.transpose()	
	A=L.multiply(L_T)

	#assertion to ensure that A and x can be multipled (Ax=b)
	assert A.columns==x.rows
	
	b=A.multiply(x)

	print "-------------------------"
	print "Expected Solution"
	print "-------------------------"
	print "x ="
	print x

	#call Choleski
	Choleski(A,b)			
				
	#At the output of a Choleski execution, A contains L and b contains y			
	# print "L = "
	# print A

	# print "y = "
	# print b	

#runs tests to ensure proper functioning of choleski fcn
#For assignment 1, Problem 1, part (b) and (c)	
def runCholeskiTests():
	L=Matrix([	[2, 0],
				[4, 3]
			])

	x=Matrix([	[3],
				[23]
			])	

	#x must be vector of length=rows of lower matrix L
	testCholeski(L,x)

	L=Matrix([	[-3, 0, 0],
				[3, 23, 0],
				[67, -89, 10]
			])

	x=Matrix([	[20],
				[-45],
				[46]
			])	

	#x must be vector of length=rows of lower matrix L
	testCholeski(L,x)

	L=Matrix([	[345.56, 0, 0, 0],
				[34, -423.5539, 0, 0],
				[2958.478, -2474.8, 10, 0],
				[123.495, 123.039, 9340.0, -2349.3]
			])

	x=Matrix([	[3746895.92783],
				[-4774975.8368],
				[29649030.764993],
				[9793479.3469802]
			])	

	#x must be vector of length=rows of lower matrix L
	testCholeski(L,x)		


#runCholeskiTests()