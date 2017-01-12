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
from lattice import Lattice
from lattice_generator import parse_yaml_and_fabric
import numpy as np
from tracks import Track

def track_design(lattice):
	vorlauf_teilchen = Proton()
	vorlauf_teilchen.out()
	dictprnt(SUMMARY,'summary track_design')

	vorlauf_spur = Track(vorlauf_teilchen)
	for pos in lattice.seq[0:28]:
		element,s0,s1 = pos
		print('{}\t(#{}, pos {:.4f}) labelled \'{}\''.format(element.__class__,id(element),s0,element.label))
		beam = element.beam
		element_matrix = element.matrix
		last = vorlauf_spur.get_last_point()
		print('track before >>',last)
		new = element_matrix.dot(last[0:6])
		last[0:6] = new
		last[6] = Beam.soll.tkin
		print('track after >>',last)
		vorlauf_spur.append_pos(s0,vorlauf_teilchen,last)
# 	vorlauf_spur.out()


def test0(filepath):
	print('\nTest0')
	lattice = parse_yaml_and_fabric(filepath)
	SUMMARY['lattice length [m]'] = CONF['lattice_length']  = lattice.length
# 	dictprnt(CONF,'CONF'); print()
	track_design(lattice)
# ---------------------------------------
if __name__ == '__main__':
	filepath = 'fodo_with_10cav_per_RF(2).yml'    ## the default input file (YAML syntax)
	test0(filepath)
