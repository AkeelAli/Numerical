import os

for i in range(3,11):
	filename = "squareCoax"+str(i)+"_D"
	os.system("python MeshSolver.py "+filename)
