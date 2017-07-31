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
import sys, os
from math import radians,sqrt,pi,degrees
import yaml

from setutil import CONF,SUMMARY,Proton,DEBUG
from setutil import objprnt,dictprnt
import elements as ELM
from lattice import Lattice

## parse and generate latttice
def lod2d(l):    ##list of dicts to dict
    return {k:v for d in l for k,v in d.items()}

def instanciate_element(item):
    # DEBUG('instanciate_element: instanciate {}'.format(item))
    key = item[0]
    attributes = item[1]
    if key == 'D':
        length   = attributes['length']
        label    = attributes['ID']
        instance =  ELM.D(length=length,label=label,particle=CONF['sollteilchen'])
    elif key == 'QF':
        length   = attributes['length']
        label    = attributes['ID']
        dBdz     = attributes["B'"]
        kq       = dBdz/CONF['sollteilchen'].brho
        instance = ELM.QF(k0=kq,length=length,label=label,particle=CONF['sollteilchen'])
    elif key == 'QD':
        length   = attributes['length']
        label    = attributes['ID']
        dBdz     = attributes["B'"]
        kq       = dBdz/CONF['sollteilchen'].brho
        instance = ELM.QD(k0=kq,length=length,label=label,particle=CONF['sollteilchen'])
    elif key == 'RFG':
        gap       = attributes['gap']
        label     = attributes['ID']
        Ez        = attributes["Ez"]
        PhiSoll   = radians(attributes["PhiSync"])
        fRF       = attributes["fRF"]
        U0        = Ez * gap
        dWf       = CONF['dWf']
        instance  =  ELM.RFG(U0=U0,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,particle=CONF['sollteilchen'],dWf=dWf)
    elif key == 'RFC':
        gap       = attributes['gap']
        length    = attributes['length']
        label     = attributes['ID']
        Ez        = attributes["Ez"]
        PhiSoll   = radians(attributes["PhiSync"])
        fRF       = attributes["fRF"]
        U0        = Ez * gap
        dWf       = CONF['dWf']
        instance  =  ELM.RFC(U0=U0,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,length=length,particle=CONF['sollteilchen'],dWf=dWf)
    elif key == 'GAP':
        gap       = attributes['gap']
        label     = attributes['ID']
        Ez        = attributes["Ez"]
        PhiSoll   = radians(attributes["PhiSync"])
        fRF       = attributes["fRF"]
        U0        = Ez * gap
        dWf       = CONF['dWf']
        instance  =  ELM.GAP(U0=U0,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,particle=CONF['sollteilchen'],dWf=dWf)
    elif key == 'MRK':
        label     = attributes['ID']
        instance  = ELM.MRK(label=label)
    else:
        raise RuntimeError('unknown element type: ',key)
    # DEBUG('instanciate_element: {} instance created'.format(label),'')
    return (label,instance)

def factory(input_file):
#--------
    def make_lattice(latticeList,segments):
        lattice = Lattice()
        # DEBUG('make_lattice for sollteilchen\n'+CONF['sollteilchen'].string())
        for segID in latticeList:        #loop segments in lattice
            # DEBUG('segID in make_lattice()',segID)
            for seg in segments:     #scan for segment in segment-definition
                if segID in seg:
                    # DEBUG('found '+segID,seg)
                    elementList = seg[segID]
                    break    #after found == true
            for element in elementList: #loop over elements in element list
                # DEBUG('element in '+segID,element)
                elementClass = element['type']
                elmItem = (elementClass,element)
                # DEBUG('elmItem in make_lattice',elmItem)
                (label,instance) = instanciate_element(elmItem)  #INSTANCIATE!!
                lattice.add_element(instance)  #add element instance to lattice
        return lattice   #the complete lattice
# --------
    def read_flags(in_data):
    #returns ==> {...}
        flags_list = in_data['flags']
        flags      = lod2d(flags_list)
        CONF['dWf'] = SUMMARY['acc. ON']                   = flags['accON']
        CONF['periodic'] = SUMMARY['ring lattice']         = flags['periodic']
        CONF['verbose']                                    = flags['verbose']
        return flags
# --------
    def read_parameters(in_data):
    #returns ==> {...}
        parameter_list = in_data['parameters']
        parameters     = lod2d(parameter_list)
        if 'frequency'        in parameters: CONF['frequenz']         = parameters['frequency']
        if 'B_grad_f'         in parameters: CONF['qf_gradient']      = parameters['B_grad_f']
        if 'B_grad_d'         in parameters: CONF['qd_gradient']      = parameters['B_grad_d']
        if 'Tkin'             in parameters: CONF['injection_energy'] = parameters['Tkin']
        if 'emitx_i'          in parameters: CONF['emitx_i']          = parameters['emitx_i']
        if 'emity_i'          in parameters: CONF['emity_i']          = parameters['emity_i']
        if 'betax_i'          in parameters: CONF['betax_i']          = parameters['betax_i']
        if 'betay_i'          in parameters: CONF['betay_i']          = parameters['betay_i']
        if 'alfax_i'          in parameters: CONF['alfax_i']          = parameters['alfax_i']
        if 'alfay_i'          in parameters: CONF['alfay_i']          = parameters['alfay_i']
        # if 'Dp/p'             in parameters: CONF['Dp/p']             = parameters['Dp/p']   # will be calculated
        if 'Ez'               in parameters: CONF['Ez_feld']          = parameters['Ez']
        if 'phi_sync'         in parameters: CONF['soll_phase']       = parameters['phi_sync']
        if 'Dz'               in parameters: CONF['Dz']               = parameters['Dz']
        if 'gap'              in parameters: CONF['spalt_laenge']     = parameters['gap']
        if 'cav_len'          in parameters: CONF['cavity_laenge']    = parameters['cav_len']
        if 'ql'               in parameters: CONF['ql']               = parameters['ql']
        if 'windings'         in parameters: CONF['n_coil']           = parameters['windings']
        CONF['wellenlänge']    = CONF['lichtgeschwindigkeit']/CONF['frequenz']
        CONF['spalt_spannung'] = CONF['Ez_feld']*CONF['spalt_laenge']
        return parameters
# --------
    def expand_reduce(in_data):
    #--------
        def read_elements(in_data):
            element_list = in_data['elements']         ## is list of dicts
            for elm in element_list:
                for elmID,attList in elm.items():      ## put key as ID in attribute dict
                    attList.append(dict(ID=elmID))
            return element_list
    #--------
        def read_segments(in_data):
            segment_list = in_data['segments']
            return segment_list
    #--------
        def read_lattice(in_data):
            lattice_segment_list= in_data['lattice']
            lattice_title = lattice_segment_list[0]['title']
            del lattice_segment_list[0]         #pull label off
            CONF['lattice_version'] = lattice_title
            return lattice_segment_list
    #--------
        def reduce_seg_def(segList):
            segs = lod2d(segList)      #{'SEG1':[...],'SEG2':[...],...}
            # DEBUG('segs in reduce_seg_def()',segs)
            segments=[]
            for segID,elmList in segs.items():
                # DEBUG(segID,elmList)
                elements=[]
                for elm in elmList:
                    elm = lod2d(elm)
                    # DEBUG('elm',elm)
                    elements.append(elm)
                segments.append({segID:elements})
            # DEBUG("segments",segments)
            return segments
        elemement_def = read_elements(in_data)
        segment_def   = read_segments(in_data)
        lattice_def   = read_lattice(in_data)
        # DEBUG('elemement_def in expand_reduce()',elemement_def)
        # DEBUG('segment_def in expand_reduce()',segment_def)
        # DEBUG('lattice_def in expand_reduce()',lattice_def)
        segments = reduce_seg_def(segment_def)
        # DEBUG('segments in expand_reduce()',segments)
        latticeList=[]
        for segSubList in lattice_def:
            nsuper = segSubList[0]
            CONF['nsuper'] = nsuper
            del segSubList[0]              #pull nsuper off
            # DEBUG('segSubList in expand_reduce()',segSubList)
            for i in range(nsuper):        #expand nsuper times
                for k in segSubList:
                    latticeList.append(k)
        return (latticeList,segments)
    ## factory body --------
    SUMMARY['input file'] = CONF['input_file'] = input_file

    with open(input_file,'r') as fileobject:
        in_data = yaml.load(fileobject)
    fileobject.close()

    read_flags(in_data)
    read_parameters(in_data)
    # DEBUG('CONF after read _parameters()',CONF.__dict__)
    (latticeList,segments) = expand_reduce(in_data)
    # DEBUG('latticeList in factory()',latticeList)      # def of all segments in lattice
    # DEBUG('segments in factory()',segments)            # def of all segments
    CONF['sollteilchen'](tkin=CONF['injection_energy'])# (WICHTIG) set sollteilchen energy
    lattice = make_lattice(latticeList,segments)
    # DEBUG('lattice_generator >> full lattice\n',lattice.string())
    SUMMARY['lattice length [m]'] = CONF['lattice_length']  = lattice.length
    # DEBUG('SUMMARY in factory()',SUMMARY)
    return lattice    #end of factory(...)

def parse_yaml_and_fabric(input_file,factory=factory):   ## delegates to factory
    return factory(input_file)

## utilities
def test0():
    print('\nTEST0')
    wfl= []
    fileobject=open('template.yml','r')
    wfl= yaml.load(fileobject)
    print(yaml.dump(wfl,default_flow_style=True))
    for i,v in iter(wfl.items()):
        print(i,' =\t',v)
    seg = wfl['segments']
    print(seg)
    print('=== segment ===')
    for i in seg:
        print(i)
    lattice = wfl['lattice']
    print('=== lattice ===')
    for l in lattice:
        for i in l:
            print(i)

def test1(input_file):
    print('\nTEST1')
    lattice = parse_yaml_and_fabric(input_file)
## main ----------
if __name__ == '__main__':
    test0()
    test1('fodo_with_10cav_per_RF(4).yml')

