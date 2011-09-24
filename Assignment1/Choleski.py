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

def Choleski(A,b):
	n=A.rows
	
	#checks to see that A is square, and b matches A
	assert n==A.columns
	assert n==b.rows
	
	#TODO checks to see if A is symmetric?
	
	#ELIMINATION
	for j in range(1,n+1):
		if (A.get(j,j)<=0): return -1
		
		A.set(j,j,math.sqrt(A.get(j,j)))
		b.set(j,1,b.get(j,1)/A.get(j,j))

		#added this to set the upper part of A (L) to 0 [optional in the algorithm]
		for i in range(1,j):
			A.set(i,j,0)
		
		for i in range(j+1,n+1):
			A.set(i,j,A.get(i,j)/A.get(j,j))
			b.set(i,1,b.get(i,1)-(A.get(i,j)*b.get(j,1)))
			
			for k in range(j+1,i+1):
				A.set(i,k,A.get(i,k)-(A.get(i,j)*A.get(k,j)))
				
	#BACK SUBSTITUTION
	x=Matrix(i=n,j=1)
	
	#counting backwards starting from n
	for i in range(n,0,-1):
		#compute the summation
		sum=0
		for j in range(i+1,n+1):
			sum+=A.get(j,i)*x.get(j,1)
		
		x.set(i,1,(b.get(i,1)-sum)/A.get(i,i))
		
	print "x = "
	print x

import math

A=Matrix([	[2,-1,0],
			[-1,2,-1],
			[0,-1,2]
		])

		
b=Matrix([	[1],
			[-1],
			[7]
		])
		

Choleski(A,b)			
			

print "L = "
print A

print "y = "
print b	

		

