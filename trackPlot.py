#!/Users/klotz/SIMULINAC_env/bin/python
# -*- coding: utf-8 -*-
__version__='v11.0.2'
"""
Copyright 2015 Wolf-Dieter Klotz <wdklotz@gmail.com>
This file is part of the SIMULINAC code

    SIMULINAC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation,Buttonnp either version 3 of the License, or
    any later version.

    SIMULINAC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SIMULINAC.  If not, see <http://www.gnu.org/licenses/>.
"""
# NOTE full with old unused or unfinished code
import numpy as NP
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.patches import Ellipse
from matplotlib.widgets import Button
from scipy.stats import norm
from math import sqrt
from setutil import DEBUG_ON,DEBUG_OFF

# NOTE below still used by bunch: from trackPlot import histPlot, poincarePlot 
def histPlot(x,mu,sigma):
    """ a historgram plot """
    import matplotlib.mlab as mlab

    num_bins = 50
    # the histogram of the data
    (n, bins, patches) = plt.hist(x, num_bins, density=1, facecolor='green', alpha=0.5)
    # add a 'best fit' line
    y = norm.pdf(bins, mu, sigma)
    plt.plot(bins, y, 'r--')
    # plt.xlabel('Smarts')
    # plt.ylabel('Probability')
    plt.title(r'$\mu$= {:.5e}, $\sigma$= {:.5e}'.format(mu,sigma))
    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)
def poincarePlot(ax,xyvalues1, xyvalues2, box, max, projections=(0,0)):
    """ 
    Scatter plot with projection histograms 
    IN:
        xyvalues    = tuple (float,float)
        box         = text 
        max         = tuple (float,float)
        projections = tuple (bool,bool)
    """
    from matplotlib.ticker import NullFormatter

    x1, y1 = xyvalues1  # sample 1
    x2, y2 = xyvalues2  # sample2
    # ax = plt.subplot(121)
    if projections == (0,0):
        # the scatter plot
        ax.scatter(x1,y1,s=1)
        if xyvalues2 != (0,0):
            ax.scatter(x2,y2,s=1,c='r')

        # set tick label size
        ax.tick_params(labelsize='xx-small')

        # place a text box in upper left in axes coords
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)  # these are matplotlib.patch.Patch properties
        ax.text(0.05, 0.95, box, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)

        (xmax,ymax) = max
        plt.xlim([-xmax,xmax])
        plt.ylim([-ymax,ymax])
        plt.autoscale(enable=False,axis='both')
        return
    
    elif projections != (0,0):
        nullfmt = NullFormatter()         # no labels
        # plt.delaxes(ax)
        # definitions for the axes
        left, width    = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h       = bottom + height + 0.02
        left_h         = left + width + 0.02

        rect_scatter = [left, bottom, width, height]
        rect_histx   = [left, bottom_h, width, 0.2]
        rect_histy   = [left_h, bottom, 0.2, height]
        
        # the scatter plot
        axScatter = plt.axes(rect_scatter)
        axScatter.scatter(x1, y1, s=1)
        if xyvalues2 != (0,0):
            axScatter.scatter(x2, y2, s=1, c='r')
        
        # set tick label size
        axScatter.tick_params(labelsize='xx-small')

        # place a text box in upper left in axes coords
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)  # these are matplotlib.patch.Patch properties
        axScatter.text(0.05, 0.95, box, transform=axScatter.transAxes, fontsize=10, verticalalignment='top', bbox=props)

        (prx,pry)   = projections
        (xmax,ymax) = max
        if prx == 1:
            # x projection
            axHistx      = plt.axes(rect_histx)
            # no labels
            axHistx.xaxis.set_major_formatter(nullfmt)
            axHistx.yaxis.set_major_formatter(nullfmt)
            # now determine nice limits by hand:
            binwidthx = xmax/100.
            limx = (int(xmax/binwidthx) + 1) * binwidthx
            axScatter.set_xlim((-limx, limx))
            binsx = NP.arange(-limx, limx + binwidthx, binwidthx)
            axHistx.hist(x1, bins=binsx, color='black')
            # axHistx.set_xlim(axScatter.get_xlim())

        if pry == 1:
            # y projection
            axHisty      = plt.axes(rect_histy)
            # no labels
            axHisty.xaxis.set_major_formatter(nullfmt)
            axHisty.yaxis.set_major_formatter(nullfmt)
            # now determine nice limits by hand:
            binwidthy = ymax/100.
            limy = (int(ymax/binwidthy) + 1) * binwidthy
            axScatter.set_ylim((-limy, limy))
            binsy = NP.arange(-limy, limy + binwidthy, binwidthy)
            axHisty.hist(y1, bins=binsy, orientation='horizontal', color='black')
            # axHisty.set_ylim(axScatter.get_ylim())

    return

# https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html
def confidence_ellipse(x, y, ax, n_std=3, facecolor='none', **kwargs):
    # """
    # Create a plot of the covariance confidence ellipse of `x` and `y`

    # Parameters
    # ----------
    # x, y : array_like, shape (n, )
    #     Input data.

    # ax : matplotlib.axes.Axes
    #     The axes object to draw the ellipse into.

    # n_std : float
    #     The number of standard deviations to determine the ellipse's radiuses.

    # Returns
    # -------
    # matplotlib.patches.Ellipse

    # Other parameters
    # ----------------
    # kwargs : `~matplotlib.patches.Patch` properties
    # """
    
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = NP.cov(x, y)
    pearson = cov[0, 1]/NP.sqrt(cov[0, 0]*cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = NP.sqrt(1 + pearson)
    ell_radius_y = NP.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),width=ell_radius_x * 2,height=ell_radius_y * 2,facecolor=facecolor,**kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = NP.sqrt(cov[0, 0]) * n_std
    mean_x = NP.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = NP.sqrt(cov[1, 1]) * n_std
    mean_y = NP.mean(y)

    transf = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)
def scatterInOut(xlive,ylive,xloss,yloss,xymax,box_txt,ax):
    ax.scatter(xlive,ylive,s=1,color='blue')
    ax.scatter(xloss,yloss,s=1,color='red')

    plt.xlim(-xymax[0],xymax[0])
    plt.ylim(-xymax[1],xymax[1])

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)
    ax.set_title(box_txt)

    # "https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html#sphx-glr-gallery-statistics-confidence-ellipse-py"
    # "https://matplotlib.org/3.1.1/gallery/statistics/confidence_ellipse.html#sphx-glr-gallery-statistics-confidence-ellipse-py"
    x=NP.concatenate((xlive,xloss))
    y=NP.concatenate((ylive,yloss))
    confidence_ellipse(x,y,ax,n_std=1,label=r"$1\sigma$",edgecolor="firebrick")
    confidence_ellipse(x,y,ax,n_std=2,label=r"$2\sigma$",edgecolor="fuchsia",linestyle="--")
    confidence_ellipse(x,y,ax,n_std=3,label=r"$3\sigma$",edgecolor="blue",linestyle=":")
    ax.legend()
    return ax

def scatter11(live,lost,abscisse,ordinate,txt):
    """
    2 scatter plots in a 11 grid
    live, lost are instances of Bunch
    abscisse, ordinate are integer coordinate indexes
    text is string
    """
    title = dict(initial=f'IN {txt}',final=f'OUT {txt}')
    loc   = dict(initial=0,final=-1)
    golden = (1.+sqrt(5.))/2.; width = 10; height = width/golden
    fig   = plt.figure(num=f'scatter plot {txt}',constrained_layout=False, figsize=(width, height))

    def plotit(*args):
        nblive   = args[0]
        nblost   = args[1]
        loc      = args[2]
        live     = args[3]
        lost     = args[4]
        abscisse = args[5]
        ordinate = args[6]
        title    = args[7]
        subplot  = args[8]
        xymax    = args[9]

        x=NP.array([]); y=NP.array([])
        nbtotal = nblive + nblost
        for particle in iter(live): # live particles
            track  = particle.track
            tpoint = track.getpoints()[loc]
            point  = tpoint()
            x = NP.append(x,point[abscisse]*1.e3)    # [mm]
            y = NP.append(y,point[ordinate]*1.e3)
        xymax0=NP.array([NP.amax(NP.abs(x)), NP.amax(NP.abs(y))])
        xlost=NP.array([]); ylost=NP.array([])
        xymax1=NP.array([0.,0.])
        if nblost != 0:     # lost particles
            for particle in iter(lost): # lost particles
                track  = particle.track
                tpoint = track.getpoints()[loc]
                point  = tpoint()
                xlost = NP.append(xlost,point[abscisse]*1.e3)    # [mm]
                ylost = NP.append(ylost,point[ordinate]*1.e3)
            xymax1=NP.array([NP.amax(NP.abs(xlost)), NP.amax(NP.abs(ylost))])
        # axis scales
        if NP.array_equal(xymax,NP.array([0.,0.])):
            xymax=NP.fmax(xymax0,xymax1)
            xymax = 1.03 * xymax   # add 3% margin
        DEBUG_OFF(xymax)

        box_text = f"{title} {nbtotal} particles"
        ax = plt.subplot(subplot)
        ax.set_title(box_text)
        plt.xlim(-xymax[0],xymax[0])
        plt.ylim(-xymax[1],xymax[1])
        ax.axvline(c='grey', lw=1)
        ax.axhline(c='grey', lw=1)

        # "https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html#sphx-glr-gallery-statistics-confidence-ellipse-py"
        # "https://matplotlib.org/3.1.1/gallery/statistics/confidence_ellipse.html#sphx-glr-gallery-statistics-confidence-ellipse-py"
        confidence_ellipse(x,y,ax,n_std=1,label=r"$1\sigma$",edgecolor="firebrick")
        confidence_ellipse(x,y,ax,n_std=2,label=r"$2\sigma$",edgecolor="fuchsia",linestyle="--")
        confidence_ellipse(x,y,ax,n_std=3,label=r"$3\sigma$",edgecolor="blue",linestyle=":")
        ax.legend()
        ax.scatter(x,y,s=1)
        if nblost !=0: 
            ax.scatter(xlost,ylost,s=1,color='red')
        return xymax
    # IN
    xymax = plotit(
        live.nbparticles(),
        lost.nbparticles(),
        loc['initial'],
        live,
        lost,
        abscisse,
        ordinate,
        title['initial'],
        121,
        NP.array([0.,0.])
        )
    # OUT
    plotit(
        live.nbparticles(),
        0,
        loc['final'],
        live,
        None,
        abscisse,
        ordinate,
        title['final'],
        122,
        xymax     # take xymax from plot before for equal axis scales
        )
    # adjust: left, bottom, right, top, wspace, hspace
    plt.subplots_adjust(wspace=0.15) 
    return

if __name__ == '__main__':
    print("what?")
