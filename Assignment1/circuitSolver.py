from Choleski import Matrix, Choleski

#makes the appropriate call to the Choleski method given 
#the incidence matrix A, matrix Y, E and J
def solveCircuit(A,Y,E,J):
	elapsedTime=Choleski(A.multiply(Y).multiply(A.transpose()),A.multiply(J.subtract(Y.multiply(E))))
	print "Time Taken: "+ str(elapsedTime)

def buildMatrices(filename):
	global Y, E, J, A
	
	#assumes proper formatting (no error detection/correction)
	f=open(filename)
	lines=f.readlines()
	f.close()

	branches=int(lines[0])
	nodes=int(lines[1])

	Y=Matrix(i=branches,j=branches)
	E=Matrix(i=branches,j=1)
	J=Matrix(i=branches,j=1)
	A=Matrix(i=nodes,j=branches)

	#NB: lines list uses array indexing (starts at 0), but Matrix class uses indexing starting at 1.

	for lineNum in range(2,(branches+1)+1): #skip first two lines already read (adjusted for array indexing, i.e. starts at 0)
		branch=lineNum-1
		for i,v in enumerate(lines[lineNum].split(',')):
			if (i==0):
				J.set(branch,1,float(v))
			elif (i==1):
				Y.set(branch,branch,1/float(v))
			else:
				E.set(branch,1,float(v))
				
	for lineNum in range((branches+1)+1,(branches+1+nodes)+1):
		node=lineNum-(branches+1)
		for j,v in enumerate(lines[lineNum].split(',')):
			A.set(node,j+1,float(v))

while (True):			
	filename=raw_input("Enter path to input file(q to quit): ")
	if (filename=='q'):
		break
	buildMatrices(filename)
	solveCircuit(A,Y,E,J)



