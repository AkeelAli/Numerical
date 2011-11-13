from Matrix import Matrix

def cg(A,b):
	#guess x all 0s
	x = Matrix(i=A.columns,j=1)

			

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

