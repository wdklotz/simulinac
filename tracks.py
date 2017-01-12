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
from setup import CONF,SUMMARY,Beam,Proton,dictprnt
# from lattice import Lattice
import numpy as np

class Track(object):
	def __init__(self,particle,particle_number=0):
		self.track_points = np.zeros(8)
		self.track_points[0] = 1.0
		self.track_points[6] = particle.tkin
		self.nbof_points = 1
		self.particle_number = particle_number

	def append_pos(self,pos,particle,new):
		track = np.zeros(8)       #new empty track point
		track += new
		track[6] = particle.tkin
		track[7] = pos
		self.track_points = np.append(self.track_points,track)
		self.nbof_points +=1
# 		print('points >>',self.track_points)

	def get_last_point(self):
		last = self.track_points.reshape(self.nbof_points,8)[-1]
		return last

	def out(self):
		points = self.track_points.reshape(self.nbof_points,8)
		for p in points:
			print('s={:.3f} tk={:.3f} x={:.3f} x\'={:.3f} y={:.3f} y\'={:.3f} z={:.3f} z\'={:.3f}'.format(p[-1],p[-2],p[0],p[1],p[2],p[3],p[4],p[5]))

