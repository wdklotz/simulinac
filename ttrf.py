import matplotlib.pyplot as plt
import numpy as np

from setutil import tblprnt

from setutil import CONF,SUMMARY,Particle,Proton,dictprnt,collect_summary,DEBUG

def test2D(X,Y,Z):
    from matplotlib.ticker import MaxNLocator
    from matplotlib.colors import BoundaryNorm
    znp = np.array(Z)
    levels = MaxNLocator(nbins=15).tick_values(znp.min(), znp.max())
    cmap = plt.get_cmap('PiYG')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    X = [X[i]*1.e3 for i in range(len(X))]
    Y = [Y[i]*1.e-6 for i in range(len(Y))]

    ## plot    
    fig = plt.figure()
    ax  = plt.subplot(111)
    pos = ax.contourf(X, Y, Z, cmap=cmap, norm=norm)
    fig.colorbar(pos, ax=ax)
    ax.set_title(title)
    ax.set_xlabel('gap [mm]')
    ax.set_ylabel('f [MHz]')
    plt.show(block=False)

def test3D(X,Y,Z):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,linewidth=0, antialiased=False)    
    fig.colorbar(surf, shrink=0.5, aspect=5)    
    ax.set_title(title)
    plt.show(block=False)
if __name__ == '__main__':
    tk = 50.
    gap = 0.044
    fRF = 816.e6
    anz = 50
    particle = Proton(tk)
    title = 'ttf for {} with tkin[MeV] {:4.1f}'.format(particle.name,particle.tkin)
    x = gapi = [0.01+i*(2.*gap/anz) for i in range(anz+1)]
    # x = [x[i]*1.e3 for i in range(len(x))]
    y = fRFi = [100.e6 + i*(fRF/anz) for i in range(anz+1)]
    X,Y = np.meshgrid(x,y)
    Z = [[particle.trtf(a,b) for a in x] for b in y]

    test2D(X,Y,Z)
    # test3D(X,Y,Z)
