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
from setup import Particle
import numpy as np
from elements import MDIM,XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO

#                    x   x'  y   y'  z   z'          Tk          1   s   1
sollStart=np.array([ 0., 0., 0., 0., 0., 0., Particle.soll.tkin, 1., 0., 1.])

class Track(object):

	#---CLASS part
	soll=None           #track of reference particle

	def string(p):   #single point to string
		s = 'x={:.3f} x\'={:.3f} y={:.3f} y\'={:.3f} z={:.3f} z\'={:.3f} s={:.3f} tk={:.5f} '.format(p[XKOO],p[XPKOO],p[YKOO],p[YPKOO],p[ZKOO],p[ZPKOO],p[SKOO],p[EKOO])
		return s

	#---INSTANCE part
	def __init__(self, particle_number=0, start=sollStart):
		self.track_points = start
		self.nbof_points = 1
		self.particle_number = particle_number

	def push(self,new):
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
			str += Track.str(p)
		return str

Track.soll = Track()      #track of Sollteilchen

