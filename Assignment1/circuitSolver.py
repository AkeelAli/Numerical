from Choleski import Matrix, Choleski

#makes the appropriate call to the Choleski method given 
#the incidence matrix A, matrix Y, E and J
def solveCircuit(A,Y,E,J):
	Choleski(A.multiply(Y).multiply(A.transpose()),A.multiply(J.subtract(Y.multiply(E))))

print "Circuit 1"
A=Matrix([	[1, -1] ])
Y=Matrix([	[0.1, 0],
			[0, 0.1]
		])		
E=Matrix([	[0],
			[10]
		])
J=Matrix([	[0],
			[0]
		])
solveCircuit(A,Y,E,J)

		
print "Circuit 2"
A=Matrix([	[1, -1] ])
Y=Matrix([	[0.1, 0],
			[0, 0.1]
		])		
E=Matrix([	[0],
			[0]
		])	
J=Matrix([	[10],
			[0]
		])
solveCircuit(A,Y,E,J)


print "Circuit 4"
A=Matrix([	[-1, 1, 1, 0],
			[0, -1, 0, 1]
		])
Y=Matrix([	[0.1, 0, 0, 0],
			[0, 0.2, 0, 0],
			[0, 0, 0.1, 0],
			[0, 0, 0, 0.2]
		])		
E=Matrix([	[10],
			[0],
			[0],
			[0]
		])	
J=Matrix([	[0],
			[0],
			[0],
			[10]
		])
solveCircuit(A,Y,E,J)

print "Circuit 5"
A=Matrix([	[-1, -1, 0, 0, 1, 0],
			[0, 1, -1, 1, 0, 0],
			[0, 0, 0, -1, -1, 1]
		])
Y=Matrix([	[0.05, 0, 0, 0, 0, 0],
			[0, 0.1, 0, 0, 0, 0],
			[0, 0, 0.033333, 0, 0, 0],
			[0, 0, 0, 0.03333, 0, 0],
			[0, 0, 0, 0, 0.1, 0],
			[0, 0, 0, 0, 0, 0.03333]
		])		
E=Matrix([	[10],
			[0],
			[0],
			[0],
			[0],
			[0]
		])	
J=Matrix([	[0],
			[0],
			[0],
			[0],
			[0],
			[0]
		])
solveCircuit(A,Y,E,J)