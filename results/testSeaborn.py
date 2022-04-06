import seaborn as sns
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def fakelog2(x):
	if x == 0:
		return 0
	else:
		return math.log2(x)

def norm(x,d):
	if(d == 0):
		return x
	else:
		return x-d

def removeZeroCol(df):
	l = []
	for x in df:
		if df[x].min() == 0 and df[x].max() == 0:
			l.append(x)
	return df.drop(l, axis=1)

def getDataFrame(filename):
	with open(filename) as fp:
	    width = len(fp.readline().strip().split(','))	
	    dtypes = {i : np.float64 for i in range(width)}
	    fp.seek(0)
	    df = pd.read_csv(fp,sep=',',header=0,dtype=dtypes)
	    df = df.dropna(axis=1)
	    df.columns = df.columns.map(np.float64)
	    return df

def plotErrBarRawAndNorm(OGdf, df):
	
	

	# #remove 0 columns
	# df = removeZeroCol(df)
	# OGdf = removeZeroCol(OGdf)

	d = {x : df[x].mean() for x in df}
	normdf = df.copy()
	for x in normdf: 
		normdf[x] = normdf[x].apply(lambda y : norm(y,OGdf[x][0]))

	# fig, axs = plt.subplots(ncols=2)
	plt.figure(figsize=(10, 6), dpi=400)
	g = sns.scatterplot(x="variable", 
						y="value", 
						data=pd.melt(df), 
						marker='_',
						# ax=axs[0],
						edgecolor=None,
						color = "blue",
						linewidth=0.5,
						label="Variants")
	# g = sns.scatterplot(data=df.T, 
	# 					markers=["_" for _ in range(len(df))],
	# 					ax=axs[0],
	# 					edgecolor=None,
	# 					linewidth=1,
	# 					palette=["blue" for _ in range(len(df))],
	# 					label="Variants")

	g = sns.scatterplot(x="variable", 
						y="value", 
						data=pd.melt(OGdf), 
						marker="_",
						# ax=axs[0],
						color = "red",
						linewidth=0.5,
						s=500,
						label="Original")
	g.set_xticks(df.columns.map(int))
	# g.set(xlabel="Weight", ylabel="log2(count)")
	g.set(xlabel="Number of Active S-boxes", ylabel="log2(count)")

	# g = sns.lineplot(x="variable", 
	# 				 y="value", 
	# 				 data=pd.melt(normdf), 
	# 				 marker=None, 
	# 				 linestyle='', 
	# 				 ci=100, 
	# 				 err_style="bars",
	# 				 err_kws={"elinewidth":2}, 
	# 				 ax=axs[1])
	# g.set_xticks(normdf.columns.map(int))
	# g.set(xlabel="Weight", ylabel="log2(countVariant) - log2(countOriginal)")

	# axs[0].legend([],[], frameon=False)

	# axs[0].set_title("Variants vs. Original")
	# axs[1].set_title("Difference of variants compared to original")
	# plt.show()

def plotStep(OGdf, df):
	plt.figure(figsize=(10, 6), dpi=400)
	g = sns.lineplot(data=df.T, 
					 marker=None,
					 linewidth=0.5,
					 label="Variants",
					 palette=["blue" for _ in range(len(df))],
					 dashes=None,
					 ci=None,
					 drawstyle='steps-post')

	# g = sns.lineplot(x="variable", 
	# 				 y="value", 
	# 				 data=pd.melt(df), 
	# 				 marker=None,
	# 				 linewidth=1,
	# 				 label="Variants",
	# 				 color="blue",
	# 				 ci=None,
	# 				 drawstyle='steps-post')
	g = sns.lineplot(x="variable", 
					 y="value", 
					 data=pd.melt(OGdf), 
					 marker=None,
					 linewidth=0.5,
					 label="Variants",
					 color="red",
					 drawstyle='steps-post')
	g.legend([],[], frameon=False)
	g.set(xlabel="Weight", ylabel="log2(count)")

filename = "data_ClusterTrail_3r_bound8_346ppiso_1iso_1Graph_noSSB"
print("Data from " + filename)

# df = getDataFrame(filename)
# df = df.applymap(fakelog2)

# OGdf = getDataFrame(filename+"_OG")
# OGdf = OGdf.applymap(fakelog2)

# plotErrBarRawAndNorm(OGdf,df)
# plt.savefig("plot_"+filename+".svg",dpi=800)

# mindf = df.min()
# maxdf = df.max()
# og = OGdf.min() #only one entry anyway

# indexlist = mindf.index
# mid = len(indexlist)//2
# if len(indexlist)%2 == 1:
# 	mid += 1
# for i in range(mid):
# 	w = indexlist[i]
# 	s = str(int(w)) + " & " + "%.2f"%mindf[w] + " & " "%.2f"%maxdf[w] + " & " + "%.2f"%og[w]
# 	if i+mid < len(indexlist):
# 		s += " & "
# 		w = indexlist[i+mid]
# 		s += str(int(w)) + " & " + "%.2f"%mindf[w] + " & " "%.2f"%maxdf[w] + " & " + "%.2f"%og[w]
# 		s += "\\\\ \\hline"
# 	else:
# 		s += "\\\\ \\cline{1-4}"
# 	print(s)

# t = [0,0,0,0]
# for i in range(len(indexlist)):
# 	w = indexlist[i]
# 	if abs(og[w] - mindf[w]) > t[0]:
# 		t = [abs(og[w] - mindf[w]), og[w], mindf[w], w]

# 	if abs(og[w] - maxdf[w]) > t[0]:
# 		t = [abs(og[w] - maxdf[w]), og[w], maxdf[w], w]

# print("Max Gap : " + str(t[0]) + " for weight " + str(t[3]) + ", " + str(t[1]) + " for OG and " + str(t[2]) + " for alt perm")

# t = [0,0,0,0]
# for i in range(len(indexlist)):
# 	w = indexlist[i]
# 	if abs(maxdf[w] - mindf[w]) > t[0]:
# 		t = [abs(maxdf[w] - mindf[w]), maxdf[w], mindf[w], w]

# print("Max Gap over all alternatives : " + str(t[0]) + " for weight " + str(t[3]) + ", " + str(t[1]) + " for max and " + str(t[2]) + " for min")

df = getDataFrame(filename)
df = df.applymap(fakelog2)

OGdf = getDataFrame(filename+"_OG")
OGdf = OGdf.applymap(fakelog2)

plotStep(OGdf,df)
plt.savefig("plot_"+filename+".svg",dpi=800)

mindf = df.min()
maxdf = df.max()
og = OGdf.min() #only one entry anyway

indexlist = mindf.index
mid = len(indexlist)//2
if len(indexlist)%2 == 1:
	mid += 1
for i in range(mid):
	w = indexlist[i]
	s = str(w) + " & " + "%.2f"%mindf[w] + " & " "%.2f"%maxdf[w] + " & " + "%.2f"%og[w]
	if i+mid < len(indexlist):
		s += " & "
		w = indexlist[i+mid]
		s += str(w) + " & " + "%.2f"%mindf[w] + " & " "%.2f"%maxdf[w] + " & " + "%.2f"%og[w]
		s += "\\\\ \\hline"
	else:
		s += "\\\\ \\cline{1-4}"
	print(s)

t = [0,0,0,0]
for i in range(len(indexlist)):
	w = indexlist[i]
	if abs(og[w] - mindf[w]) > t[0]:
		t = [abs(og[w] - mindf[w]), og[w], mindf[w], w]

	if abs(og[w] - maxdf[w]) > t[0]:
		t = [abs(og[w] - maxdf[w]), og[w], maxdf[w], w]

print("Max Gap : " + str(t[0]) + " for weight " + str(t[3]) + ", " + str(t[1]) + " for OG and " + str(t[2]) + " for alt perm")

t = [0,0,0,0]
for i in range(len(indexlist)):
	w = indexlist[i]
	if abs(maxdf[w] - mindf[w]) > t[0]:
		t = [abs(maxdf[w] - mindf[w]), maxdf[w], mindf[w], w]

print("Max Gap over all alternatives : " + str(t[0]) + " for weight " + str(t[3]) + ", " + str(t[1]) + " for max and " + str(t[2]) + " for min")