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

from setutil import Particle,DEBUG
from elements import MDIM,XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
import elements as ELM

class Track(object):    #is an ordered list of track-points. A track-point is an array of MDIM coordinates.

	def __init__(self, particle_number=0, start=None):
		self.track_points = start
		self.particle_number = particle_number
		self.nbof_points = 1

	def nbPoints(self):
		return self.nbof_points

	def append(self,new):
		self.track_points = np.append(self.track_points,new)
		self.nbof_points +=1

	def points(self):
		return self.track_points.reshape(self.nbof_points,MDIM)

	def point_at(self,n):
		return self.points()[n]

	def first(self):
		first = self.point_at(0)
		return first

	def last(self):
		last = self.point_at(-1)
		return last

	def first_str(self):
		return Track.string(self.first())

	def last_str(self):
		return Track.string(self.last())

	def points_str(self):
		points = self.points()
		str = ''
		for p in points:
			str += Track.string(p)
		return str

	def string(p):   #single point to string
		s = 'x={:.3e} x\'={:.3e} y={:.3e} y\'={:.3e} z={:.3e} z\'={:.3e}  tk={:.5f} s={:.3f} '.format(p[XKOO],p[XPKOO],p[YKOO],p[YPKOO],p[ZKOO],p[ZPKOO],p[EKOO],p[SKOO])
		return s
	soll = None

#default track-point    x   x'  y   y'  z   z'          Tk          1   s   1
Track.soll = Track(start=np.array([ 0., 0., 0., 0., 0., 0., Particle.soll.tkin, 1., 0., 1.]))

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
	"""
	Call in a loop to create terminal progress bar
	@params:
		iteration   - Required  : current iteration (Int)
		total       - Required  : total iterations (Int)
		prefix      - Optional  : prefix string (Str)
		suffix      - Optional  : suffix string (Str)
		decimals    - Optional  : positive number of decimals in percent complete (Int)
		length      - Optional  : character length of bar (Int)
	"""
	percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
	filledLength = int(length * iteration // total)
	bar = fill * filledLength + '-' * (length - filledLength)
	print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
	# Print New Line on Complete
	if iteration == total:
		print()
#
# Sample Usage of printProgressBar(...)
#
#
# from time import sleep
#
# # make a list
# items = list(range(0, 57))
# i = 0
# l = len(items)
#
# # Initial call to print 0% progress
# printProgressBar(i, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
# for item in items:
#     # Do stuff...
#     sleep(0.1)
#     # Update Progress Bar
#     i += 1
#     printProgressBar(i, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
#
# # Sample Output
# Progress: |█████████████████████████████████████████████-----| 90.0% Complete

def trackSoll(lattice):
	"""
	Tracks the reference particle through the lattice and redefines the lattice element parameters to
	adapt them to the energy of the accellerated reference particle.
	"""
	soll_track = Track.soll       #track of reference particle
	for ipos in lattice.seq:
		element,s0,s1 = ipos
# 		DEBUG('\n{}\t(#{}, pos {:.4f}) label \'{}\''.format(element.__class__,id(element),s0,element.label))
		ti = soll_track.last()                #track: at entrance
# 		DEBUG('\t\ti >>',Track.string(ti))
		element.adapt_for_energy(ti[EKOO])    #enery adopt
		tf = element.matrix.dot(ti)           #track: at exit
		soll_track.append(tf)                 #append
# 		deltaE = tf[EKOO] - ti[EKOO]
# 		DEBUG('\t\tf >>',Track.string(tf),' deltaE[KeV] >>',deltaE*1.e3)
# 	DEBUG('complete track\n{}'.format(soll_track.points_string()))
	DEBUG('{}'.format(soll_track.first_str()))
	DEBUG('{}'.format(soll_track.last_str()))

def track(lattice,bunch):
	"""
	Tracks a bunch of particles through the lattice
	lattice: a list of elements a.k.a. _matrix'es
	bunch: a list of independent Tracks
	"""
	from time import sleep
	printProgressBar(0, bunch.nbTracks(), prefix = 'Progress:', suffix = 'Complete', length = 50)
	for (count,particle_track) in enumerate(bunch.tracks()):
		ti = particle_track.first()
		for ipos in lattice.seq:
			element,s0,s1 = ipos
			tf = element.matrix.dot(ti)      #track through!
			if isinstance(element,ELM.MRK):
				particle_track.append(tf)
# 				deltaE = tf[EKOO] - ti[EKOO]
# 				DEBUG('\t\tf >>',Track.string(tf),' deltaE[KeV] >>',deltaE*1.e3)
			ti = tf
		sleep(1.0e-3)
		printProgressBar(count, bunch.nbTracks(), prefix = 'Progress:', suffix = 'Complete', length = 50)
# 		DEBUG('complete track\n{}'.format(particle_track.points_str()))
# 		DEBUG('FIRST: {}'.format(particle_track.first_str()))
# 		DEBUG('{} LAST: {}'.format(count,particle_track.last_str()))

