reset()

def getStructSSB(M):

	inputDone = set()
	structSSB = []
	while len(inputDone) < 16:
		for i in range(16):
			if i not in inputDone:
				inset = set([i])
				outset = set()
				prev_inset = set()
				prev_outset = set()

				while(prev_inset != inset or prev_outset != outset):
					prev_inset = set(inset)
					prev_outset = set(outset)
					for j in inset:
						for k in range(16):
							if M[j][k] == 1:
								outset.add(k)

					for j in outset:
						for k in range(16):
							if M[k][j] == 1:
								inset.add(k)

				structSSB.append(len(inset))
				for x in inset:
					inputDone.add(x)
				break

	return tuple(sorted(structSSB))

def genGraphFile(M, filename):
	f = open(filename, "w")
	for i in range(M.nrows()):
		for j in range(M.ncols()):
			if M[i][j] != 0:
				f.write(str(i)+">"+str(j)+"\n")
	f.close()

load("allCentralDiGraphs.sage") #AllDiGraphs list as integer matrices

listDiGraphsMatrix = dict()

for Mdata in AllDiGraphs:
	M = Matrix(16,16)
	for i in range(16):
		rowData = Mdata[i]
		for j in range(16):
			if(rowData & (1 << j) != 0):
				M[i,j] = 1
	structSSB = getStructSSB(M)
	if structSSB in listDiGraphsMatrix:
		listDiGraphsMatrix[structSSB].append(M)
	else:
		listDiGraphsMatrix[structSSB] = [M]

# save(listDiGraphsMatrix, "listDiGraphsMatrix")

#DiGraphs that lead to no SSB
L = load("listDiGraphsMatrix.sobj")[tuple([16])]
L = listDiGraphsMatrix[tuple([16])]

dirName = 'graphFiles'
if not os.path.exists(dirName):
    os.mkdir(dirName)
    print("Directory " , dirName ,  " Created ")
else:    
    print("Directory " , dirName ,  " already exists")

for i in range(len(L)):
	genGraphFile(L[i], "./graphFiles/graphNoSSB"+str(i))

L = listDiGraphsMatrix[tuple([4,4,4,4])]
for i in range(len(L)):
	genGraphFile(L[i], "./graphFiles/graph4444SSB"+str(i))
