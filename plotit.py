import pandas as pd
import matplotlib.pyplot as plt
plt.style.use("ggplot")

sir_out = pd.read_csv("sir_out",sep=" ",header=None,names=["t","S","E","I","A","H","R","D"],index_col=False)

sline = plt.plot("t","S","",data=sir_out,color="red",linewidth=2)
eline = plt.plot("t","E","",data=sir_out,color="green",linewidth=2)
iline = plt.plot("t","I","",data=sir_out,color="purple",linewidth=2)
aline = plt.plot("t","A","",data=sir_out,color="pink",linewidth=2)
hline = plt.plot("t","H","",data=sir_out,color="yellow",linewidth=2)
rline = plt.plot("t","R","",data=sir_out,color="blue",linewidth=2)
dline = plt.plot("t","D","",data=sir_out,color="black",linewidth=2)
plt.xlabel("Time",fontweight="bold")
plt.ylabel("Number",fontweight="bold")
legend = plt.legend(title="Population",loc=10,bbox_to_anchor=(1.25,0.5))
frame = legend.get_frame()
frame.set_facecolor("white")
frame.set_linewidth(0)

plt.show()