
import subprocess
import sys

if __name__ == "__main__":

	out = subprocess.check_output(["./glasgow-subgraph-solver/glasgow_subgraph_solver", "partialGraph", "./graphFiles/"+sys.argv[1], "--print-all-solutions"])
	out = out.decode('utf-8')

	Liso = []
	while out.find("mapping") != -1:
		out = out[out.find("mapping"):]
		row = out[:out.find("\n")]
		out = out.replace(row+"\n", "")

		row = row.replace("mapping = ", "")
		row = row.replace(") (", ",")
		row = row.replace(")", "")
		row = row.replace("(", "")
		row = row.split(",")
		row = [[int(y) for y in x.split(" -> ")] for x in row]
		
		L = [None for _ in range(16)]
		for x,y in row:
			L[x] = y
		buff = ""
		for x in L:
			if x == None:
				buff += "-1"
			else:
				buff += str(x)
			buff += ","
		print(buff.rstrip(","))