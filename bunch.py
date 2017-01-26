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
from math import sqrt
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from setutil import DEBUG,CONF,Particle
from elements import MDIM,XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
from tracks import Track

def sigmas(alfa,beta,epsi):
	gamma = (1.+ alfa**2)/beta
	sigma  = sqrt(epsi*beta)
	sigmap = sqrt(epsi*gamma)
	return sigma,sigmap

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

def poincarePlot(x,y,whazit,sctrplt,max=(0.,0.)):       #scatter plot

# 	max values
	xmax = max[0]
	if xmax == 0.:
		xmax = np.max(np.fabs(x))
	ymax = max[1]
	if ymax == 0.:
		ymax = np.max(np.fabs(y))
# 	DEBUG('xmax,ymax in poincarePlot() >>',xmax,ymax)

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

	if xmax != 0:
		# now determine nice binning limits by hand
		binwidthx=xmax/100.
		limx = (int(xmax/binwidthx)+1)*binwidthx
		binsx = np.arange(-limx, limx + binwidthx, binwidthx)
		axHistx.hist(x, bins=binsx,normed=1)
		axHistx.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')

	if ymax != 0.:
		# now determine nice binning limits by hand
		binwidthy=ymax/100.
		limy = (int(ymax/binwidthy)+1)*binwidthy
		binsy = np.arange(-limy, limy + binwidthy, binwidthy)
		axHisty.hist(y, bins=binsy, normed=1, orientation='horizontal')
		axHisty.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')

	return xmax,ymax

class EmittanceContour(object):
	def twiss_conjugate(x,alfa,beta,epsi):
		gamma = (1.+ alfa**2)/beta
		a=beta
		b=2.*alfa*x
		c=gamma*x**2-epsi
		d=(b**2-4.*a*c)
		d=sqrt(d)
		xp1=(-b+d)/(2.*a)
		xp2=(-b-d)/(2.*a)
		return (xp1,xp2)

	def __init__(self,nbof_tracks,args):
		sigx,sigxp = sigmas(CONF['alfax_i'],CONF['betax_i'],CONF['emitx_i'])
		if args['random']:
			Xrand = sigx*(2.*np.random.random_sample((nbof_tracks,))-1.)
		else:
			Xrand = np.linspace(-sigx*(1.-1.e-3),sigx*(1.-1.e-3),nbof_tracks)
		X=[]; XP=[]
		for x in Xrand:
			xp1,xp2 = EmittanceContour.twiss_conjugate(x,CONF['alfax_i'],CONF['betax_i'],CONF['emitx_i'])
			X.append(x)
			XP.append(xp1)
			X.append(x)
			XP.append(xp2)
		sigy,sigyp = sigmas(CONF['alfay_i'],CONF['betay_i'],CONF['emity_i'])
		if args['random']:
			Yrand = sigy*(2.*np.random.random_sample((nbof_tracks,))-1.)
		else:
			Yrand = np.linspace(-sigy+1.e-5,sigy-1.e-5,nbof_tracks)
		Y=[]; YP=[]
		for y in Yrand:
			yp1,yp2 = EmittanceContour.twiss_conjugate(y,CONF['alfay_i'],CONF['betay_i'],CONF['emity_i'])
			Y.append(y)
			YP.append(yp1)
			Y.append(y)
			YP.append(yp2)
		tk_in = Particle.soll.tkin                           #energy at entrance
# 		DEBUG('X >>', X)
# 		DEBUG('XP >>',XP)
# 		DEBUG('Y >>', Y)
# 		DEBUG('YP >>',YP)
		self.tracklist=[]           #all Tracks in a bunch
		for i in range(2*nbof_tracks):
			start=np.array([ 0., 0., 0., 0., 0., 0., tk_in, 1., 0., 1.])
			start[XKOO]=X[i]
			start[XPKOO]=XP[i]
			start[YKOO]=Y[i]
			start[YPKOO]=YP[i]
			self.tracklist.append(Track(particle_number=i,start=start))
# 			DEBUG(self.tracklist[-1].first_str())
# 			DEBUG(self.tracklist[-1].last_str())

class Gauss1D(object):
	def __init__(self,nbof_tracks,args):
		sigx,sigxp = sigmas(CONF['alfax_i'],CONF['betax_i'],CONF['emitx_i'])
		sigy,sigyp = sigmas(CONF['alfay_i'],CONF['betay_i'],CONF['emity_i'])
		sigz  = CONF['Dz']
		sigzp = CONF['Dp/p']
		X  = sigx  *np.random.randn(nbof_tracks)     #gauss distribution X
		XP = sigxp *np.random.randn(nbof_tracks)
		Y  = sigy  *np.random.randn(nbof_tracks)
		YP = sigyp *np.random.randn(nbof_tracks)
		Z  = sigz *np.random.randn(nbof_tracks)
		ZP = sigzp *np.random.randn(nbof_tracks)
		plane = args['plane']
		tk_in = Particle.soll.tkin                           #energy at entrance
# 		DEBUG('X >>', X)
# 		DEBUG('XP >>',XP)
# 		DEBUG('Y >>', Y)
# 		DEBUG('YP >>',YP)
# 		DEBUG('Z  >>',Z)
# 		DEBUG('ZP >>',ZP)
		self.tracklist=[]           #all Tracks in a bunch
		for i in range(nbof_tracks):
			start=np.array([ 0., 0., 0., 0., 0., 0., tk_in, 1., 0., 1.])
			if plane[0]:
				start[XKOO]=X[i]
			if plane[1]:
				start[XPKOO]=XP[i]
			if plane[2]:
				start[YKOO]=Y[i]
			if plane[3]:
				start[YPKOO]=YP[i]
			if plane[4]:
				start[ZKOO]=Z[i]
			if plane[5]:
				start[ZPKOO]=ZP[i]
			self.tracklist.append(Track(particle_number=i,start=start))
# 			DEBUG(self.tracklist[-1].first_str())
# 			DEBUG(self.tracklist[-1].last_str())

class Bunch(object):  #is a list of Tracks, which is a list of track-points, which is an array of coordinates
	def __init__(self,init=True):
		self.nbof_tracks = 750
		self.tracklist = []     #all Tracks in a bunch -- a list []
		self.plane = (1,1,1,1,1,1)
		self.distclass=Gauss1D
		if init: self.initPhaseSpace({'plane':self.plane})
	#---
	def nbTracks(self):        #nbof tracks per bunch
		return self.nbof_tracks
	#---
	def set_nbTracks(self,nb):
		self.nbof_tracks = nb
	#---
	def tracks(self):          #all tracks
		return self.tracklist
	#---
	def last(self):
		return self.tracklist[-1]
	#---
	def nbPointsPTrack(self):     #nbof points per track
		return self.tracklist[-1].nbPoints()
	#---
	def set_plane(self,plane):
		self.plane = plane
	#---
	def set_distClass(self,distclass):
		self.distclass = distclass
	#---
	def initPhaseSpace(self,args):
		self.tracklist = self.distclass(self.nbof_tracks,args).tracklist
		if self.distclass == EmittanceContour:
			self.set_nbOfParticles(len(self.tracklist))  #EmittanceCountour doubles nbof-tracks

def test1(alfa,beta,epsi):
	N = 20000
	sigma = sqrt(epsi*beta)
	gamma = (1.+ alfa**2)/beta
	sigmap = sqrt(epsi*gamma)
	x  = sigma *np.random.randn(N)
	xp = sigmap *np.random.randn(N)
	xp = xp*xp/(sigmap*sigmap)

	DEBUG('x >>',x)
	DEBUG('x\' >>',xp)

# 	plt.figure()
# 	h1 = plt.subplot2grid((2,1),(0,0))
# 	histPlot(x,0.,sigma)
# 	h2 = plt.subplot2grid((2,1),(1,0))
# 	histPlot(xp,0.,sigmap)

	plt.figure()
	sp = plt.subplot()
	poincarePlot(x,xp,'x-x\'',sp)
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
