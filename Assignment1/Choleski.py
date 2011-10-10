import math
import copy
from time import clock

class MatrixElement:
	
	def __init__(self,row,column,value):
		self.value=value
		self.i=row
		self.j=column
	
	def __repr__(self):
		return "("+str(self.i)+","+str(self.j)+")="+str(self.v)
	
	
class Matrix:

	#constructor can be called 2 ways
	#1) sending a list of lists of values
	#2) sending i and j to build an empty matrix of size i,j
	def __init__(self,listOfLists=None,i=None,j=None):
		self.elements={}
		self.rows=None
		self.columns=None
	
		if (listOfLists is None and (not i is None) and (not j is None)):
			self.rows=i
			self.columns=j
			
			for j in range(1,self.columns+1):
				for i in range(1,self.rows+1):
					self.elements[str(i)+","+str(j)]=MatrixElement(i,j,0)
			
		elif (not listOfLists is None):
			#TODO assert that listOfLists is a list
			#assert isinstance(listOfLists, List)
			
			self.rows=len(listOfLists)
			assert self.rows>0 #there needs to be at least one list inside the listOfLists
			
			#check that all the lists have same number of items & assign it to self.columns
			for i in range(0,self.rows):
				if (self.columns==None): #first iteration
					self.columns=len(listOfLists[i])
				else:
					assert len(listOfLists[i])==self.columns
				
			#create the matrix elements
			for i,list in enumerate(listOfLists):
				for j,value in enumerate(list):
					self.elements[str(i+1)+","+str(j+1)]=MatrixElement(i+1,j+1,value)
	
	#returns Matrix product of self with passed multiplier Matrix
	def multiply(self,multiplier):
		#return None if inner dimensions don't match
		if (self.columns!=multiplier.rows):
			return None
		
		innerDimension=self.columns
		
		result=Matrix(i=self.rows,j=multiplier.columns)
		
		for j in range(1, multiplier.columns+1):
			for i in range(1, self.rows+1):
				sum=0
				for k in range(1, innerDimension+1):
					sum+=self.get(i,k)*multiplier.get(k,j)
				
				result.set(i,j,sum)

		return result
	
	#returns the self's transpose as a Matrix
	def transpose(self):
		result=Matrix(i=self.columns,j=self.rows)
		
		for j in range(1, self.columns+1):
			for i in range(1, self.rows+1):
				result.set(j,i,self.get(i,j))
				
		return result
	
	def subtract(self,subtrahend):
		if (self.rows!=subtrahend.rows or self.columns!=subtrahend.columns):
			return None
		
		result=Matrix(i=self.rows,j=self.columns)
		
		for j in range(1, self.columns+1):
			for i in range(1, self.rows+1):
				result.set(i,j,self.get(i,j)-subtrahend.get(i,j))
				
		return result

	def get(self,i,j):
		return self.elements[str(i)+","+str(j)].value
		
	def set(self,i,j,value):
		self.elements[str(i)+","+str(j)].value=value
	
	#sets all values to 0
	def clear(self):
		for j in range(1,columns+1):
			for i in range(1, rows+1):
				self.set(i,j,0)
	
	def __repr__(self):
		string=""
		#using a list for the values that will be added to the string to be 
		#printed according to our preferred format (i.e. %.2f)
		string_values=[]
		for i in range(1,self.rows+1):
			for j in range(1,self.columns+1):
				string+=" "
				
				#negative numbers have one extra dash character (-)
				#less one decimal place for -ve numbers to keep alignment
				value=self.elements[str(i)+","+str(j)].value
				if (value>=0):
					string+="%.3f"
				else:
					string+="%.2f"
				
				string_values.append(value)
			string+="\n"
		return string % tuple(string_values)

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
	print "-------------------------"
	print "Solution"
	print "-------------------------"
	print "x ="
	print x
	# print "=============================================="
	# print "END Choleski Function Output						"
	# print "=============================================="
	print ""
	
	return elapsedTime

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