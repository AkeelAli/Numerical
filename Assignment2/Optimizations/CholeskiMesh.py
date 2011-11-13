from Choleski import Choleski
from Matrix import Matrix
import MeshStructure
import pickle
import sys

def printResults(nodeH,pList,filename):
	f = open(filename,"w")
		
	f.write("node,x,y,potential,predefined\n")
	for i in range(1,len(nodeH)+1):
		f.write(str(nodeH[i].number)+","+str(nodeH[i].x)+","+str(nodeH[i].y)+","+str(nodeH[i].voltage))
		if ( nodeH[i].number in pList):
			f.write(",y\n")
		else:
			f.write("\n")
		
	f.close()

def getAlpha(x1,y1,x2,y2,x3,y3,x,y):
	denominator = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)
	if (denominator == 0):
		return 99999 #for infinity
	else:
		return ((y2-y3)*(x-x3)+(x3-x2)*(y-y3)) / (denominator*1.0)
	
def getBeta(x1,y1,x2,y2,x3,y3,x,y):
	denominator = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)
	if (denominator == 0):
		return 99999 #for infinity
	else:	
		return ((y3-y1)*(x-x3)+(x1-x3)*(y-y3)) / (denominator*1.0)
	
filename = sys.argv[1]
print "--------------------------------"
print "Creating results csv for File "+filename+".msh"
print "--------------------------------"
#load data
print "#load data"
fd = open(filename+".dump",'r')
data = pickle.load(fd)
fd.close()
A = data[0]
b = data[1]
nodeH = data[2]
pList = data[3]
U_con_order = data[4]

#Call Choleski
print "#Call Choleski"
U_f = Choleski(A,b)

#put results in appropriate nodes in nodeH
for i in range(1,U_f.rows+1):
	nodeH[U_con_order[i-1]].voltage = U_f.get(i,1)

#nodeH at this point has all the proper voltages

#for each triangle, see if the point we're after is inside
MeshStructure.buildStructs(filename+".msh") #we need triangleL
triangleL=MeshStructure.triangleL
for i in range(0,len(triangleL)):
	n1 = triangleL[i].node1
	n2 = triangleL[i].node2
	n3 = triangleL[i].node3
	alpha = getAlpha(n1.x,n1.y,n2.x,n2.y,n3.x,n3.y,0.05,0.05)
	beta = getBeta(n1.x,n1.y,n2.x,n2.y,n3.x,n3.y,0.05,0.05)
	gamma = 1 - alpha - beta
	
	if ((alpha >= 0 and alpha <=1) and (beta >=0 and beta <=1) and (gamma >=0 and gamma <=1)):
		potential = nodeH[n1.number].voltage*alpha + nodeH[n2.number].voltage*beta + nodeH[n3.number].voltage*gamma
		print "free nodes = "+str(U_f.rows)+" , (0.05,0.05) = "+str(potential)
	
	
#print results
#print "#print results"
#printResults(nodeH,pList,filename+".csv")

