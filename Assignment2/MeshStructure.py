#3 data structures for Problem 2 (a)
nodeH = {}
triangleH = {}
pList = []


class Node:
	def __init__(self,number,x,y,voltage):
		self.x = x
		self.y = y
		self.number = number
		self.voltage = voltage
		
	def __repr__(self):
		return str(self.number)+"("+str(self.x)+","+str(self.y)+")="+str(self.voltage)
		
class Triangle:
	def __init__(self,n1,n2,n3):
		self.node1 = n1
		self.node2 = n2
		self.node3 = n3
		
	def area(self):
		pass
		
	def __repr__(self):
		return str(self.node1)+"|"+str(self.node2)+"|"+str(self.node3)
		

def buildStructs(filename):
	import re
	
	f = open(filename,'r')
	
	# Create Node Hash#
	global nodeH
	numNodes = int(f.readline())
	
	for i in range(1,numNodes + 1):
		line = f.readline()
		m = re.match('(\d+.\d+) (\d+.\d+)',line)
		
		nodeH[i] = Node(i,float(m.group(1)), float(m.group(2)), 0.0)
	
	
	# Create Triangle Hash#
	global triangleH
	numT = int(f.readline())
	
	for i in range(1, numT + 1):
		line = f.readline()
		#TODO what is the last value?
		m = re.match('(\d+) (\d+) (\d+) (\d+.\d+)',line)
		nodeNum1 = int(m.group(1))
		nodeNum2 = int(m.group(2))
		nodeNum3 = int(m.group(3))

		#Store each triangle 6 times in the hash (3*2*1) for ease of access:(x,yz) or (x,z,y) or (y,x,z) etc.
		for j in (str(nodeNum1),str(nodeNum2),str(nodeNum3)):
			for k in (str(nodeNum1),str(nodeNum2),str(nodeNum3)):
				for w in (str(nodeNum1),str(nodeNum2),str(nodeNum3)):
					if (not (j == k or j == w or k == w)):
						triangleH["("+j+","+k+","+w+")"] = \
							Triangle(nodeH[nodeNum1],nodeH[nodeNum2],nodeH[nodeNum3])
	

	# Read in boundary conditions #
	# Update Nodes + store prescribed in list #
	global pList
	numBC = int(f.readline())
	
	for i in range(1,numBC + 1):
		line = f.readline()
		
		m = re.match('(\d+) (\d+.\d+)',line)
		nodeNum = int(m.group(1))
		voltage = float(m.group(2))
		
		nodeH[nodeNum].voltage = voltage;
		pList.append(nodeNum)
	
	
	f.close()
	