from Choleski import Choleski
from Matrix import Matrix
import MeshStructure
import pickle
import sys

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

#print results
print "#print results"
printResults(nodeH,pList,filename+".csv")
