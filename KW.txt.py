import numpy as np
import matplotlib.pyplot as plt
from math import radians

x=np.linspace(0,10,100)
print(x)
fig=plt.figure()
ax=plt.axes()
ax.plot(x,np.cos(x)+np.sin(x+radians(30.)),'--')
fig.show()
