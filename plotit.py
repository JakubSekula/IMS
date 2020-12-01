import pandas as pd
import matplotlib.pyplot as plt
plt.style.use("ggplot")

sir_out = pd.read_csv("sir_out",sep=" ",header=None,names=["t","S","I","D","A","R","T","H"],index_col=False)

#sline = plt.plot("t","S","",data=sir_out,color="red",linewidth=1)
iline = plt.plot("t","S","",data=sir_out,color="blue",linewidth=1)
dline = plt.plot("t","I","",data=sir_out,color="red",linewidth=1)
aline = plt.plot("t","D","",data=sir_out,color="green",linewidth=1)
rline = plt.plot("t","A","",data=sir_out,color="black",linewidth=1)
tline = plt.plot("t","R","",data=sir_out,color="blue",linewidth=1,linestyle='dashed')
hline = plt.plot("t","T","",data=sir_out,color="red",linewidth=1,linestyle='dashed')
eline = plt.plot("t","H","",data=sir_out,color="green",linewidth=1,linestyle='dashed')

plt.xlabel("Time",fontweight="bold")
plt.ylabel("Number",fontweight="bold")
legend = plt.legend(title="Population",loc=10,bbox_to_anchor=(1.25,0.5))
frame = legend.get_frame()
frame.set_facecolor("white")
frame.set_linewidth(0)

plt.show()