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
from xml_utils.XmlDataAdaptor import XmlDataAdaptor
from setutil import FLAGS,PARAMS,SUMMARY
from lattice_generator import factory
from lattice import Lattice
print(dir())

def makeSecDa(root_da,sec):
    pass

input_file = '25_09_2017_versuche_70_200MeV'
input_dir  = 'yml/'
input_ext  = '.yml'
input      = '{}{}{}'.format(input_dir,input_file,input_ext)
lattice = factory(input)

root_da  = XmlDataAdaptor(name='Alceli')
sections = lattice.get_sections()       #sections is a [Section,...]

for section in sections:
    if len(section.seq) == 0: continue
    sec     = section.get_name()
    sec_da  = root_da.createChild(sec)
    cavs_da = sec_da.createChild('Cavities')
    for node in section.seq:
        node_da = sec_da.createChild('accElement')
        node_da.setValue('length',node.length)
        node_da.setValue('pos',node.position[1])
        node_da.setValue('name',node.label)
        node_da.setValue('type',node.__class__.__name__)
print(root_da.makeXmlText())
root_da.writeToFile('{}{}'.format(input_file,'.xml'))
