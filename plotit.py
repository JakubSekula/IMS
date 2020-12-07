import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

plt.style.use("ggplot")

sir_out = pd.read_csv("ims_out",sep=" ",header=None,names=["t","D"],index_col=False)

#iline = plt.plot("t","S","",data=sir_out,color="blue",linewidth=1, label='Cumulative infected' )
#dline = plt.plot("t","I","",data=sir_out,color="red",linewidth=1, label='Current total infected' )
#aline = plt.plot("t","DE","",data=sir_out,color="green",linewidth=1, label='Recovered' )
rline = plt.plot("t","D","",data=sir_out,color="black",linewidth=1, label='Deaths' )
#tline = plt.plot("t","R","",data=sir_out,color="blue",linewidth=1,linestyle='dashed', label='Diagnosed cumulative infected' )
#hline = plt.plot("t","T","",data=sir_out,color="red",linewidth=1,linestyle='dashed', label='Diagnosed current total infected' )
#eline = plt.plot("t","H","",data=sir_out,color="green",linewidth=1,linestyle='dashed', label='Diagnosed recovered' )


plt.xlabel("Time (days)", fontsize="13" )
plt.ylabel("Deaths", fontsize="13")

plt.legend(loc="upper left", fontsize="10")

x1,x2,y1,y2 = plt.axis()

#plt.axis((0,350,0,0.015))


plt.show()