import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from math import pi,degrees,sin,cos

columns = ["x","xp","y","yp","z","zp","T",1,"S","1"]

n = 40000
p=[]
index=[]
for i in range(n+1):
	arg = 2*pi*i/n
	x=sin(arg)
	xp=cos(arg)
	y=sin(2*arg)
	yp=cos(2*arg)
	p.append([x,xp,y,yp,0.,0.,0.,1.,0.,1.])
	index.append(arg)

df=pd.DataFrame(p,columns=columns,index=index)
df.loc[:,["x","xp","y","yp"]].plot()
pl.show()


