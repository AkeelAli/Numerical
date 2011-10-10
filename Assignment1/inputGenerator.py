from Choleski import Matrix

branches=[]

n=int(raw_input("Enter N: "))
output=raw_input("Enter output filename: ")

def populateBranches():
	global branches
	
	
	upperNodes=[1]
	node=1 #will hold the node number
	
	#build a list of branches (where a branch is identified by (node1,node2))
	#method used (start from top, going down level by level)
	for i in range(2,(n+1)+1):
		currentNodes=[]
		
		for j in range(1,i+1):
			node+=1
			currentNodes.append(node)
		
		#branches linking upperNodes level to currentNodes level
		for index,upperNode in enumerate(upperNodes):
			branches.append((upperNode,currentNodes[index]))
			branches.append((upperNode,currentNodes[index+1]))
		
		#branches linking currentNodes level
		for index,currentNode in enumerate(currentNodes):
			if (index+1==len(currentNodes)): #stop at last node
				break
			branches.append((currentNode,currentNodes[index+1]))
			
		upperNodes=currentNodes




numBranches=(3*n**2+3*n)/2
numNodes=(n**2+3*n+2)/2
populateBranches()

numBranches+=1 #add 1 to the number of branches (for source branch)

#build incidence matrix from branches list
A=Matrix(i=numNodes,j=numBranches) 
for j,branch in enumerate(branches):
	A.set(branch[0],j+1,1)
	A.set(branch[1],j+1,-1)

#set source branch
A.set(1,numBranches,-1)
A.set(numNodes,numBranches,1)

	
f=open(output,'w')

f.write(str(numBranches)+'\n')
#number of nodes (subtract one taken as ground)
f.write(str(numNodes-1)+'\n')


#print J,R,E for each branch (0,1,0)
for i in range(1,numBranches):
	f.write('0,1,0\n')
#print source branch
f.write('0,1,1\n') #1V source, 1ohm R_test
	
#print incidence matrix	
for i in range(1,A.rows): #stop at rows, and not rows+1 as we discard last node taken as ground
	f.write(str(int(A.get(i,1))))
	for j in range(2,A.columns+1):
		f.write(','+str(int(A.get(i,j))))
	f.write('\n')
	
f.close()

