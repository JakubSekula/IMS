import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

plt.style.use("ggplot")

sir_out = pd.read_csv("ims_out",sep=" ",header=None,names=["t","S","I","D","A","R","T","H"],index_col=False)

#sline = plt.plot("t","S","",data=sir_out,color="red",linewidth=1)
iline = plt.plot("t","S","",data=sir_out,color="blue",linewidth=1, label='Infected ND AS',)
dline = plt.plot("t","I","",data=sir_out,color="red",linewidth=1, label='Infected D AS',  )
aline = plt.plot("t","D","",data=sir_out,color="green",linewidth=1, label='Infected ND S' )
rline = plt.plot("t","A","",data=sir_out,color="purple",linewidth=1, label='Infected D S' )
tline = plt.plot("t","R","",data=sir_out,color="brown",linewidth=1,label='Infected D IC' )


plt.xlabel("Time (days)", fontsize="13" )
plt.ylabel("Cases (fraction of population)", fontsize="13" )

plt.legend(loc="upper right", fontsize="10")

x1,x2,y1,y2 = plt.axis()

plt.axis((0,350,0,0.0011))

plt.show()