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
from setup import CONF,SUMMARY,dictprnt
from lattice import Lattice
from lattice_generator import parse_yaml_and_fabric
import numpy as np
from tracks import Track

def track_design(lattice):
	soll_spur = Track.soll
# 	for ipos in lattice.seq[0:28]:
	for ipos in lattice.seq[0:-1]:
		element,s0,s1 = ipos
# 		print('{}\t(#{}, pos {:.4f}) labelled \'{}\''.format(element.__class__,id(element),s0,element.label))
		element_matrix = element.matrix
		before = soll_spur.last_out()
# 		print('i >>',Track.out(before))
		after = element_matrix.dot(before)
# 		print('f >>',Track.out(after))
		soll_spur.push(s0,after)
	print(soll_spur.all_out())


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
