import matplotlib.pyplot as plt
import math as m
import numpy as np

twopi = 2.*m.pi

x      = np.linspace(0,3*twopi,num=200)

y      = np.array([m.sin(a) for a in x])
rms_y  = np.array([m.sqrt(a*a) for a in y])

py      = np.array([m.sin(a+m.degrees(15.)) for a in x])
rms_py  = np.array([m.sqrt(a*a) for a in py])


# fig = ax = plt.figure()
fig,((ax11,ax12,ax13),(ax21,ax22,ax23)) = plt.subplots(nrows=2,ncols=3)
ax11.plot(x,y)
ax21.plot(x,rms_y)
ax12.plot(x,py)
ax22.plot(x,rms_py)

ax13.plot(y,py)
ax23.plot(rms_y,py)

plt.show()