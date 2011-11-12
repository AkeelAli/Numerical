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
					string+="%.5f"
				else:
					string+="%.5f"
				
				string_values.append(value)
			string+="\n"
		return string % tuple(string_values)