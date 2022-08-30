#!/Users/klotz/SIMULINAC_env/bin/python
# -*- coding: utf-8 -*-
"""
Copyright 2015 Wolf-Dieter Klotz <wdklotz@gmail.com>
This file is part of the SIMULINAC code

    SIMULINAC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    SIMULINAC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SIMULINAC.  If not, see <http://www.gnu.org/licenses/>.
"""
#TODO full with old unused or unfinished code
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import pprint, inspect

def PRINT_PRETTY(obj):
    file = inspect.stack()[0].filename
    print('DEBUG_ON ==============>  '+file)
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

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

def scatter_hist(x,y, ax, ax_histx,ax_histy):
   # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)

    # now determine nice limits by hand:
    binwidth = 0.25
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation='horizontal')

def poincare_hist(fig,ax,x,y):
    # definitions for the axes
    left, width    = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing        = 0.005


    rect_scatter = [left, bottom, width, height]
    rect_histx   = [left, bottom + height + spacing, width, 0.2]
    rect_histy   = [left + width + spacing, bottom, 0.2, height]

    # start with a square Figure
    # fig = plt.figure(figsize=(8, 8))

    # Add a gridspec with two rows and two columns and a ratio of 2 to 7 between
    # the size of the marginal axes and the main axes in both directions.
    # Also adjust the subplot parameters for a square plot.
    gs = fig.add_gridspec(2, 2,  width_ratios=(7, 2), height_ratios=(2, 7), left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.05, hspace=0.05)

    ax1       = ax.add_subplot(gs[1, 0])
    ax_histx = ax.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = ax.add_subplot(gs[1, 1], sharey=ax)
    # ax = fig.add_axes(rect_scatter)
    # ax_histx = fig.add_axes(rect_histx, sharex=ax)
    # ax_histy = fig.add_axes(rect_histy, sharey=ax)

    # use the previously defined function
    scatter_hist(x, y, ax1, ax_histx, ax_histy)

#TODO: generalize for many xyvalues
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
            binsx = np.arange(-limx, limx + binwidthx, binwidthx)
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

if __name__ == '__main__':
    print("baaaaaaa - nothing todo")
