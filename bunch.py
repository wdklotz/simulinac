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
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from setup import DEBUG,CONF,Particle
from math import sqrt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from elements import MDIM,XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
from tracks import Track

def histPlot(x,mu,sigma):      #histogram
	num_bins = 50
	# the histogram of the data
	(n, bins, patches) = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
	# add a 'best fit' line
	y = mlab.normpdf(bins, mu, sigma)
	plt.plot(bins, y, 'r--')
# 	plt.xlabel('Smarts')
# 	plt.ylabel('Probability')
	plt.title(r'$\mu$= {:.5e}, $\sigma$= {:.5e}'.format(mu,sigma))
	# Tweak spacing to prevent clipping of ylabel
	plt.subplots_adjust(left=0.15)

def sectionPlot(x,y,whazit,sctrplt):       #scatter plot

# 	max values
	xmax = np.max(np.fabs(x))
	ymax = np.max(np.fabs(y))

	# the scatter plot
	sctrplt.scatter(x,y,s=1)

	# set tick label size
	sctrplt.tick_params(labelsize='xx-small')

	# place a text box in upper left in axes coords
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)  # these are matplotlib.patch.Patch properties
	sctrplt.text(0.05, 0.95, whazit, transform=sctrplt.transAxes, fontsize=10, verticalalignment='top', bbox=props)

	# 	create new axes on the right and on the top
	divider = make_axes_locatable(sctrplt)
	axHistx = divider.append_axes('top',   1.2, pad=0.2, sharex=sctrplt)
	axHisty = divider.append_axes('right', 1.2, pad=0.2, sharey=sctrplt)
	axHistx.tick_params(labelsize='xx-small')
	axHisty.tick_params(labelsize='xx-small')

	# make some labels invisible
	plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(), visible=False)

	# now determine nice binning limits by hand
	binwidthx=xmax/100.
	binwidthy=ymax/100.
	limx = (int(xmax/binwidthx)+1)*binwidthx
	limy = (int(ymax/binwidthy)+1)*binwidthy
	binsx = np.arange(-limx, limx + binwidthx, binwidthx)
	binsy = np.arange(-limy, limy + binwidthy, binwidthy)

	# do the histograms
	axHistx.hist(x, bins=binsx)
	axHisty.hist(y, bins=binsy, orientation='horizontal')

#	axHistx.axis['bottom'].major_ticklabels.set_visible(False)
	for tl in axHistx.get_xticklabels():
		tl.set_visible(False)
	axHistx.set_yticks([0,50,100])

#	axHisty.axis['left'].major_ticklabels.set_visible(False)
	for tl in axHisty.get_yticklabels():
		tl.set_visible(False)
	axHisty.set_xticks([0,50,100])

class Bunch(object):  #is a list of Tracks, which is a list of track-points, which is an array of coordinates
	def __init__(self,nbpart=1000,init=True):
		self.nbof_particles = nbpart
		self.tracklist = None
		if init: self.initPhaseSpace()
	#---
	def nb_particles(self):
		return self.nbof_particles
	#---
	def tracks(self):
		return self.tracklist
	#---
	def initPhaseSpace(self,plane=(1,1,1,1,0,0)):
		def sigmas(alfa,beta,epsi):
			gamma = (1.+ alfa**2)/beta
			sigma  = sqrt(epsi*beta)
			sigmap = sqrt(epsi*gamma)
			return sigma,sigmap
		#---
		sigx,sigxp = sigmas(CONF['alfax_i'],CONF['betax_i'],CONF['emitx_i'])
		sigy,sigyp = sigmas(CONF['alfay_i'],CONF['betay_i'],CONF['emity_i'])
		X  = sigx  *np.random.randn(self.nbof_particles)     #gauss distribution X
		XP = sigxp *np.random.randn(self.nbof_particles)
		Y  = sigy  *np.random.randn(self.nbof_particles)
		YP = sigyp *np.random.randn(self.nbof_particles)
		tk_in = Particle.soll.tkin                           #energy at entrance
# 		DEBUG('X >>', X)
# 		DEBUG('XP >>',XP)
# 		DEBUG('Y >>', Y)
# 		DEBUG('YP >>',YP)
		self.tracklist=[]           #all Tracks in a bunch
		for i in range(self.nbof_particles):
			start=np.array([ 0., 0., 0., 0., 0., 0., tk_in, 1., 0., 1.])
			if plane[0]:
				start[XKOO]=X[i]
			if plane[1]:
				start[XPKOO]=XP[i]
			if plane[2]:
				start[YKOO]=Y[i]
			if plane[3]:
				start[YPKOO]=YP[i]
			self.tracklist.append(Track(particle_number=i,start=start))
# 			DEBUG(self.tracklist[-1].first_str())
# 			DEBUG(self.tracklist[-1].last_str())
		return self.tracklist

def test2(nbpart):
	bunch = Bunch(nbpart)
	bunch.initPhaseSpace((1,0,0,0,0,0))

def test1(alfa,beta,epsi):
	N = 20000
	sigma = sqrt(epsi*beta)
	gamma = (1.+ alfa**2)/beta
	sigmap = sqrt(epsi*gamma)
	x  = sigma *np.random.randn(N)
	xp = sigmap *np.random.randn(N)

# 	DEBUG('x >>',x)
# 	DEBUG('x\' >>',xp)

# 	plt.figure()
# 	h1 = plt.subplot2grid((2,1),(0,0))
# 	histPlot(x,0.,sigma)
# 	h2 = plt.subplot2grid((2,1),(1,0))
# 	histPlot(xp,0.,sigmap)

	plt.figure()
	sp = plt.subplot()
	sectionPlot(x,xp,'x-x\'',sp)
	plt.show(block=True)

def test0(mu,sigma):
	# example data
# 	mu = 0  # mean of distribution
# 	sigma = 1  # standard deviation of distribution
	x = mu + sigma * np.random.randn(2000)
	histPlot(x,mu,sigma)
	plt.show(True)

if __name__ == '__main__':
# 	test0(2.,1.)
	test1(CONF['alfax_i'],CONF['betax_i'],CONF['emitx_i'])
# 	test2(10)
