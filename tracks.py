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
from setup import MDIM,Particle
import numpy as np

sollStart=np.array([0.,0.,0.,0.,0.,0.,Particle.soll.tkin,1.,0.,1.])
sollStart=np.array([1.,0.,1.,0.,0.,0.,Particle.soll.tkin,1.,0.,1.])

class Track(object):
	#---CLASS part
	soll=None           #track of reference particle
	def out(p):   #single point out
		str = 's={:.3f} tk={:.3f} x={:.3f} x\'={:.3f} y={:.3f} y\'={:.3f} z={:.3f} z\'={:.3f}'.format(p[-2],p[-4],p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[8])
		return str
	#---INSTANCE part
	def __init__(self, particle_number=0, start=sollStart):
		self.track_points = start
		self.nbof_points = 1
		self.particle_number = particle_number

	def push(self,pos,new):
		self.track_points = np.append(self.track_points,new)
		self.nbof_points +=1
# 		print('points >>',self.track_points)

	def last_out(self):
		last = self.track_points.reshape(self.nbof_points,MDIM)[-1]
		return last

	def all_out(self):
		points = self.track_points.reshape(self.nbof_points,MDIM)
		str = ''
		for p in points:
			str += 's={:.3f} tk={:.3f} x={:.3f} x\'={:.3f} y={:.3f} y\'={:.3f} z={:.3f} z\'={:.3f}\n'.format(p[-2],p[-4],p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[8])
		return str

Track.soll = Track()      #track of Sollteilchen

