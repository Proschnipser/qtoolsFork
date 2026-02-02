from nexus import matrix2nexus
import numpy as np
matrix=np.matrix([[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]])
b = [matrix[x, :x + 1] for x in range(3)]

matrix2nexus(matrix,taxa=["mac", "human","rat","giraffe"],nexusfile="./test.nex",imgformat='png',plot_now=False,cmdfile=False,cmdmode='w',splitstree_location='~/splitstree4/SplitsTree')