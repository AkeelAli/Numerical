from Matrix import Matrix
import sys
import pickle
import copy
import math

LIMITERROR = 0.0000001

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

def cg(A,b):
	#guess x all 0s
	x = Matrix(i=A.columns,j=1)

	#set r and p
	r = b.subtract(A.multiply(x))
	p = b.subtract(A.multiply(x))
		
	while(math.sqrt(r.transpose().multiply(r).get(1,1)) > LIMITERROR):
		p_t = p.transpose()
		alpha = p_t.multiply(r).get(1,1) /  p_t.multiply(A.multiply(p)).get(1,1)
		
		
		# a_p for alpha*p
		a_p = copy.deepcopy(p)
		for i in range(1,a_p.rows+1):
			for j in range(1,a_p.columns+1):
				a_p.set(i,j,a_p.get(i,j)*alpha)
		
		x = x.add(a_p)
		
		r = b.subtract(A.multiply(x))
		beta = -1 * (p_t.multiply(A.multiply(r)).get(1,1) / p_t.multiply(A.multiply(p)).get(1,1))

		
		# b_p for beta*p
		b_p = p #no need to copy (as we update cell by cell of p to b_p and we don't need it later)
		for i in range(1,b_p.rows+1):
			for j in range(1,b_p.columns+1):
				b_p.set(i,j,p.get(i,j)*beta)
		p = r.add(b_p)

		
	return x
	
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

U_f = cg(A,b)

printResults(nodeH,pList,filename+"2.csv")

