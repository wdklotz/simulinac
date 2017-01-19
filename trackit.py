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
import matplotlib.pyplot as plt
import time

from lattice_generator import parse_yaml_and_fabric
from bunch import Bunch,sectionPlot
from tracks import Track,track,trackSoll
from elements import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
from setup import DEBUG

def scatterplot(bnch,xko,yko,txt):
	x=[]; y=[]
	for t in bnch.tracks():
		x.append(t.last()[xko])
		y.append(t.last()[yko])
	plt.figure()
	sp = plt.subplot()
	sectionPlot(x,y,'{} {} particles'.format(txt,bnch.nb_particles()),sp)

def test0():
	from math import radians
	from elements import RFC,RFG
	from setup import Particle
	import numpy as np
	from elements import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO

	print('\ntest0: correctness of RFG and RFC')
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

	tracks = Track(start=Track.soll.first())
	count = 10
	while count != 0:
		ti = tracks.last()
		rf4 = RFC(U0=U0 ,PhiSoll=PhiSoll ,fRF=fRF, label='rf4', particle=Particle.soll, gap=gap, length=gap, dWf=1)
		rf4.adapt_for_energy(ti[EKOO])
		tf = rf4.matrix.dot(ti)
		tracks.append(tf)
		deltaW= tf[EKOO] - ti[EKOO]
		print('track through rf4')
		print(rf4.string())
		print('ti >>',ti)
		print('tf >>',tf,' deltaW[KeV] >>',deltaW*1.e3)
		count -= 1

def test1(filepath):
	from setup import CONF,dictprnt
	print('\ntest1: trackSoll(...)')
	lattice = parse_yaml_and_fabric(filepath)
# 	SUMMARY['lattice length [m]'] = CONF['lattice_length']  = lattice.length
	dictprnt(CONF,'CONF'); print()
	trackSoll(lattice)

def trackit(filepath):
	print('\ntrackit')
	t0 = time.clock()
	lattice = parse_yaml_and_fabric(filepath)
	t1 = time.clock()
	bunch = Bunch(500)
	t2 = time.clock()
	trackSoll(lattice)
	t3 = time.clock()
	track(lattice,bunch)
	t4 = time.clock()
	scatterplot(bunch,XKOO,XPKOO,'x-x\'')
	scatterplot(bunch,YKOO,YPKOO,'y-y\'')
	scatterplot(bunch,XKOO,YKOO,'x-y')
	t5 = time.clock()
	DEBUG('total time     >> {:6.3f} [sec]'.format(t5-t0))
	DEBUG('parse lattice  >> {:6.3f} [sec] {:4.1f} [%]'.format((t1-t0),(t1-t0)/(t5-t0)*1.e2))
	DEBUG('generate bunch >> {:6.3f} [sec] {:4.1f} [%]'.format((t2-t1),(t2-t1)/(t5-t0)*1.e2))
	DEBUG('track design   >> {:6.3f} [sec] {:4.1f} [%]'.format((t3-t2),(t3-t2)/(t5-t0)*1.e2))
	DEBUG('track bunch    >> {:6.3f} [sec] {:4.1f} [%]'.format((t4-t3),(t4-t3)/(t5-t0)*1.e2))
	DEBUG('fill plots     >> {:6.3f} [sec] {:4.1f} [%]'.format((t5-t4),(t5-t4)/(t5-t0)*1.e2))
	plt.show(block=True)

# ---------------------------------------
if __name__ == '__main__':
	filepath = 'fodo_with_10cav_per_RF(2).yml'    ## the default input file (YAML syntax)
# 	test0()
# 	test1(filepath)
	trackit(filepath)
