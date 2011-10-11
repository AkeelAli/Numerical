outLength=0.2
inWidth=0.08
inHeight=0.04
ERROR=10**-5

outVoltage=0
inVoltage=10

#default h and w
h=0.02
w=1

#maximum coordinate of i (where i starts at 0)
maxI=(outLength/2)/h

maxJ=(outLength/2)/h

#maximum pts when j=maxJ where we are still on neumann boundary (and not V=10)
maxI_Neumann=(outLength/2-inWidth/2)/h   -1 # -1 to remove the boundary point (which is also 10V)
maxJ_Neumann=(outLength/2-inHeight/2)/h    -1

#will contain potential values (addressable by coordinate points)
hash={}

def buildHash():
	global hash
	
	for i in range(0,int(maxI+1)):
		for j in range(0,int(maxJ+1)):
			if (i==0 or j==0):
				hash[(i,j)]=float(outVoltage)
			# elif (i==maxI and j<=maxJ_Neumann):
				# hash[(i,j)]=0
			# elif (j==maxJ and i<=maxI_Neumann):
				# hash[(i,j)]=0
			elif (i>maxI_Neumann and j>maxJ_Neumann):
				hash[(i,j)]=float(inVoltage)
			else:
				hash[(i,j)]=float(0)

def printHash():
	import sys

	for j in range(int(maxJ),-1,-1):
		for i in range(0, int(maxI)+1):
			sys.stdout.write(str(hash[(i,j)]))
			sys.stdout.write(',')
		print ''

def solveFDM():
	iterations=0
	repeat=True
	
	while(repeat):
		iterations+=1
	
		#working bottom->up, left->right to respect order of computations as per the algorithm
		#skip i=0 , j=0 as these are the outer boundary (dirichlet)
		for i in range(1,int(maxI)+1):
			for j in range (1, int(maxJ)+1):
				if (i>maxI_Neumann and j>maxJ_Neumann):
					continue #skip the inner conductor with fixed potential 10V
				elif (i==maxI):
					hash[(i,j)]=( (1-w) * hash[(i,j)] ) + float(w)/4 * (  hash[(i-1,j)] + hash[(i,j-1)] + hash[(i-1,j)] + hash[(i,j+1)]  )#neumann
				elif (j==maxJ):
					hash[(i,j)]=( (1-w) * hash[(i,j)] ) + float(w)/4 * (  hash[(i-1,j)] + hash[(i,j-1)] + hash[(i+1,j)] + hash[(i,j-1)]  ) #neumann
				else:
					hash[(i,j)]= ( (1-w) * hash[(i,j)] ) + float(w)/4 * (  hash[(i-1,j)] + hash[(i,j-1)] + hash[(i+1,j)] + hash[(i,j+1)]  )

		
		
		#check if residual for each node is small enough
		repeat=False
		
		for i in range(1,int(maxI)+1):
			for j in range (1, int(maxJ)+1):
				if (i>maxI_Neumann and j>maxJ_Neumann):
					continue #skip the inner conductor with fixed potential 10V
				elif (i==maxI):
					continue #TODO skip neumann residual checking?
				elif (j==maxJ):
					continue #TODO skip neumann residual checking?
				else:
					residual = (  hash[(i-1,j)] + hash[(i,j-1)] + hash[(i+1,j)] + hash[(i,j+1)]  ) - 4 * hash[(i,j)]
					
				if (residual > ERROR):
					repeat = True
			
			if (repeat):
				break
		
	return iterations

def Jacobi():
	import copy
	iterations=0
	repeat=True
	
	while(repeat):
		iterations+=1
		
		hashPrevious=copy.deepcopy(hash)
		
		#working bottom->up, left->right to respect order of computations as per the algorithm
		#skip i=0 , j=0 as these are the outer boundary (dirichlet)
		for i in range(1,int(maxI)+1):
			for j in range (1, int(maxJ)+1):
				if (i>maxI_Neumann and j>maxJ_Neumann):
					continue #skip the inner conductor with fixed potential 10V
				elif (i==maxI):
					hash[(i,j)]= 0.25 * (  hashPrevious[(i-1,j)] + hashPrevious[(i,j-1)] + hashPrevious[(i-1,j)] + hashPrevious[(i,j+1)]  )#neumann
				elif (j==maxJ):
					hash[(i,j)]= 0.25 * (  hashPrevious[(i-1,j)] + hashPrevious[(i,j-1)] + hashPrevious[(i+1,j)] + hashPrevious[(i,j-1)]  ) #neumann
				else:
					hash[(i,j)]= 0.25 * (  hashPrevious[(i-1,j)] + hashPrevious[(i,j-1)] + hashPrevious[(i+1,j)] + hashPrevious[(i,j+1)]  )

		
		
		#check if residual for each node is small enough
		repeat=False
		
		for i in range(1,int(maxI)+1):
			for j in range (1, int(maxJ)+1):
				if (i>maxI_Neumann and j>maxJ_Neumann):
					continue #skip the inner conductor with fixed potential 10V
				elif (i==maxI):
					continue #TODO skip neumann residual checking?
				elif (j==maxJ):
					continue #TODO skip neumann residual checking?
				else:
					residual = (  hash[(i-1,j)] + hash[(i,j-1)] + hash[(i+1,j)] + hash[(i,j+1)]  ) - 4 * hash[(i,j)]
					
				if (residual > ERROR):
					repeat = True
			
			if (repeat):
				break
		
	return iterations
	
#Question 3(b)
# for i in range(10,20,1):
	# w=float(i)/10
	# print w
	# hash={}
	# buildHash()
	# print solveFDM()
	# print hash[(0.06/h,0.04/h)]
	# print ''
	
#Question 3(c)
# w=1.3

# for i in range(0,5):
	# h=0.02/(2**i)
	# print h
	
	# #maximum coordinate of i (where i starts at 0)
	# maxI=(outLength/2)/h

	# maxJ=(outLength/2)/h

	# #maximum pts when j=maxJ where we are still on neumann boundary (and not V=10)
	# maxI_Neumann=(outLength/2-inWidth/2)/h   -1 # -1 to remove the boundary point (which is also 10V)
	# maxJ_Neumann=(outLength/2-inHeight/2)/h    -1
	
	# hash={}
	# buildHash()
	# print solveFDM()
	# print hash[(0.06/h,0.04/h)]
	# print ''
	
#Question 3(d)
for i in range(0,5):
	h=0.02/(2**i)
	print h
	
	#maximum coordinate of i (where i starts at 0)
	maxI=(outLength/2)/h

	maxJ=(outLength/2)/h

	#maximum pts when j=maxJ where we are still on neumann boundary (and not V=10)
	maxI_Neumann=(outLength/2-inWidth/2)/h   -1 # -1 to remove the boundary point (which is also 10V)
	maxJ_Neumann=(outLength/2-inHeight/2)/h    -1
	
	hash={}
	buildHash()
	print Jacobi()
	print hash[(0.06/h,0.04/h)]
	print ''
	