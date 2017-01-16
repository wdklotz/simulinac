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
from math import radians
from setup import CONF,SUMMARY,dictprnt,Particle,DEBUG
from lattice import Lattice
from lattice_generator import parse_yaml_and_fabric
import numpy as np
from tracks import Track
from elements import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
from elements import RFC,RFG

def track_soll(lattice):
	soll_track = Track.soll       #track of reference particle
# 	for ipos in lattice.seq[0:3*28]:
	for ipos in lattice.seq:
		element,s0,s1 = ipos
# 		DEBUG('\n{}\t(#{}, pos {:.4f}) label \'{}\''.format(element.__class__,id(element),s0,element.label))
		ti = soll_track.last()
# 		DEBUG('\t\ti >>',Track.string(ti))
		element.adapt_for_energy(ti[EKOO])
		tf = element.matrix.dot(ti)      #track through!
		soll_track.push(tf)
		deltaE = tf[EKOO] - ti[EKOO]
# 		DEBUG('\t\tf >>',Track.string(tf),' deltaE[KeV] >>',deltaE*1.e3)
# 	DEBUG('complete track\n{}'.format(soll_track.points_string()))
# 	DEBUG('{}'.format(soll_track.first_str()))
# 	DEBUG('{}'.format(soll_track.last_str()))

def test0():
	print('\nTest0: correctness of RFG and RFC')
	U0=0.05; PhiSoll=radians(-20.); fRF=816.e6; gap=0.05; length=gap; dWf=1
			#RFG
# 	rf  = RFG(U0=U0 ,PhiSoll=PhiSoll ,fRF=fRF, label='rf', particle=Particle.soll, gap=gap, dWf=1)
# 	rf1 = RFG(U0=U0 ,PhiSoll=PhiSoll ,fRF=fRF, label='rf1', particle=Particle(tkin=200.), gap=gap, dWf=1)
# 	rf.string(0)
# 	rf1.string(0)
			#RFC
	rf2 = RFC(U0=U0 ,PhiSoll=PhiSoll ,fRF=fRF, label='rf2', particle=Particle.soll, gap=gap, length=gap, dWf=1)
	rf3 = RFC(U0=U0 ,PhiSoll=PhiSoll ,fRF=fRF, label='rf3', particle=Particle(tkin=200.), gap=gap, length=gap, dWf=1)
	print(rf2.string())
	print(rf3.string())

	rf2.adapt_for_energy(tkin=200.)
	print('rf2.adapt_for_energy(tkin=200.)')
	print(rf2.string())
	m = rf3.matrix - rf2.matrix
	print('rf3.matrix - rf2.matrix. Must be all zeros!\n',m)


	trackStart=np.array([ 0., 0., 0., 0., 0., 0., Particle(tkin=100.).tkin, 1., 0., 1.])
	tracks = Track(start=trackStart)
	count = 10
	while count != 0:
		ti = tracks.last()
		rf4 = RFC(U0=U0 ,PhiSoll=PhiSoll ,fRF=fRF, label='rf4', particle=Particle.soll, gap=gap, length=gap, dWf=1)
		rf4.adapt_for_energy(ti[EKOO])
		tf = rf4.matrix.dot(ti)
		tracks.push(tf)
		deltaW= tf[EKOO] - ti[EKOO]
		print('track through rf4')
		print(rf4.string())
		print('ti >>',ti)
		print('tf >>',tf,' deltaW[KeV] >>',deltaW*1.e3)
		count -= 1



def test1(filepath):
	print('\nTest1')
	lattice = parse_yaml_and_fabric(filepath)
# 	SUMMARY['lattice length [m]'] = CONF['lattice_length']  = lattice.length
# 	dictprnt(CONF,'CONF'); print()
	track_soll(lattice)
# ---------------------------------------
if __name__ == '__main__':
	filepath = 'fodo_with_10cav_per_RF(2).yml'    ## the default input file (YAML syntax)
	test0()
# 	test1(filepath)
