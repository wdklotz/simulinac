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
import numpy as np
import matplotlib.pyplot as plt

# DEBUGGING
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass

def histPlot(x,mu,sigma):
    """ a historgram plot """
    import matplotlib.mlab as mlab

    num_bins = 50
    # the histogram of the data
    (n, bins, patches) = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
    # add a 'best fit' line
    y = mlab.normpdf(bins, mu, sigma)
    plt.plot(bins, y, 'r--')
    # plt.xlabel('Smarts')
    # plt.ylabel('Probability')
    plt.title(r'$\mu$= {:.5e}, $\sigma$= {:.5e}'.format(mu,sigma))
    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)

def poincarePlot(good, whazit, max, bad=([],[]), projections=(0,0)):
    """ a scatter plot with projection histograms """
    from matplotlib.ticker import NullFormatter

    x, y = good
    ax = plt.subplot()
    if projections == (0,0):
        # the scatter plot
        # ax.scatter(x,y,s=40,color=['b','g','r','c','m','y'])
        ax.scatter(x,y,s=1)

        # set tick label size
        ax.tick_params(labelsize='xx-small')

        # place a text box in upper left in axes coords
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)  # these are matplotlib.patch.Patch properties
        ax.text(0.05, 0.95, whazit, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)

        (xmax,ymax) = max
        plt.xlim([-xmax,xmax])
        plt.ylim([-ymax,ymax])
        plt.autoscale(enable=False,axis='both')
    
    elif projections != (0,0):
        nullfmt = NullFormatter()         # no labels
        plt.delaxes(ax)
        # definitions for the axes
        left, width    = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h       = bottom + height + 0.02
        left_h         = left + width + 0.02

        rect_scatter = [left, bottom, width, height]
        rect_histx   = [left, bottom_h, width, 0.2]
        rect_histy   = [left_h, bottom, 0.2, height]
        
        xi, yi = bad              # bad are lost tracks
        # the scatter plot:
        axScatter = plt.axes(rect_scatter)
        axScatter.scatter(x, y, s=1)
        axScatter.scatter(xi, yi, s=1, color='r')

        # set tick label size
        axScatter.tick_params(labelsize='xx-small')

        # place a text box in upper left in axes coords
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)  # these are matplotlib.patch.Patch properties
        axScatter.text(0.05, 0.95, whazit, transform=axScatter.transAxes, fontsize=10, verticalalignment='top', bbox=props)

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
            axHistx.hist(x, bins=binsx, color='black')
            axHistx.set_xlim(axScatter.get_xlim())

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
            binsy = np.arange(-limy, limy + binwidthy, binwidthy)
            axHisty.hist(y, bins=binsy, orientation='horizontal', color='black')
            axHisty.set_ylim(axScatter.get_ylim())

    return

if __name__ == '__main__':
    print("baaaaaaa - nothing todo")