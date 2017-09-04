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

from setutil import PARAMS,FLAGS,SUMMARY,Proton,DEBUG,objprnt,dictprnt,zellipse
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
        instance =  ELM.D(length=length,label=label,particle=PARAMS['sollteilchen'])
    elif key == 'QF':
        length   = attributes['length']
        label    = attributes['ID']
        dBdz     = attributes["B'"]
        kq       = dBdz/PARAMS['sollteilchen'].brho
        instance = ELM.QF(k0=kq,length=length,label=label,particle=PARAMS['sollteilchen'])
    elif key == 'QD':
        length   = attributes['length']
        label    = attributes['ID']
        dBdz     = attributes["B'"]
        kq       = dBdz/PARAMS['sollteilchen'].brho
        instance = ELM.QD(k0=kq,length=length,label=label,particle=PARAMS['sollteilchen'])
    elif key == 'RFG':
        gap       = attributes['gap']
        label     = attributes['ID']
        Ez        = attributes["Ez"]
        PhiSoll   = radians(attributes["PhiSync"])
        fRF       = attributes["fRF"]
        U0        = Ez * gap
        dWf       = FLAGS['dWf']
        instance  =  ELM.RFG(U0=U0,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,particle=PARAMS['sollteilchen'],dWf=dWf)
    elif key == 'RFC':
        gap       = attributes['gap']
        length    = attributes['length']
        label     = attributes['ID']
        Ez        = attributes["Ez"]
        PhiSoll   = radians(attributes["PhiSync"])
        fRF       = attributes["fRF"]
        U0        = Ez * gap
        dWf       = FLAGS['dWf']
        instance  =  ELM.RFC(U0=U0,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,length=length,particle=PARAMS['sollteilchen'],dWf=dWf)
    elif key == 'GAP':
        gap       = attributes['gap']
        label     = attributes['ID']
        Ez        = attributes["Ez"]
        PhiSoll   = radians(attributes["PhiSync"])
        fRF       = attributes["fRF"]
        U0        = Ez * gap
        dWf       = FLAGS['dWf']
        instance  =  ELM.GAP(U0=U0,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,particle=PARAMS['sollteilchen'],dWf=dWf)
    elif key == 'MRK':
        label     = attributes['ID']
        instance  = ELM.MRK(label=label)
    else:
        raise RuntimeError('unknown element type: ',key)
    # DEBUG('instanciate_element: {} instance created'.format(label),'')
        sys.exit(1)
    try:     ## sections are not mandatory
        instance.set_section(sec=attributes['sec'])
    except:
        pass
        # instance.set_section(sec='undef')
    return (label,instance)

def factory(input_file):
#--------
    def make_lattice(latticeList,segments):
        lattice = Lattice()
        # DEBUG('make_lattice for sollteilchen\n'+PARAMS['sollteilchen'].string())
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
        flags = lod2d(flags_list) if flags_list != None else {}
        if 'accON' in flags and flags['accON']: 
            FLAGS['dWf']     = 1.
            SUMMARY['accON'] = True
        if 'periodic'    in flags: FLAGS['periodic'] = SUMMARY['ring lattice']     = flags['periodic']
        if 'egf'         in flags: FLAGS['egf']      = SUMMARY['emittance growth'] = flags['egf']
        if 'sigma'       in flags: FLAGS['sigma']    = SUMMARY['sigma tracking']   = flags['sigma']
        if 'map'         in flags: FLAGS['map']      = SUMMARY['track with map']   = flags['map']
        if 'KVprint'     in flags: FLAGS['KVprint']                                = flags['KVprint']
        if 'verbose'     in flags: FLAGS['verbose']                                = flags['verbose']
        return flags
# --------
    def read_sections(in_data):
    #returns ==> [...]
        try:     ## sections are not mandatory
            sec_list = in_data['sections'][0]
        except:
            sec_list = []
        PARAMS['sections'] = sec_list
        return sec_list
# --------
    def read_parameters(in_data):
    #returns ==> {...}
        parameter_list = in_data['parameters']
        parameters     = lod2d(parameter_list)
        if 'frequency'        in parameters: PARAMS['frequenz']         = parameters['frequency']
        if 'B_grad_f'         in parameters: PARAMS['qf_gradient']      = parameters['B_grad_f']
        if 'B_grad_d'         in parameters: PARAMS['qd_gradient']      = parameters['B_grad_d']
        if 'Tkin'             in parameters: PARAMS['injection_energy'] = parameters['Tkin']
        if 'emitx_i'          in parameters: PARAMS['emitx_i']          = parameters['emitx_i']
        if 'emity_i'          in parameters: PARAMS['emity_i']          = parameters['emity_i']
        if 'betax_i'          in parameters: PARAMS['betax_i']          = parameters['betax_i']
        if 'betay_i'          in parameters: PARAMS['betay_i']          = parameters['betay_i']
        if 'alfax_i'          in parameters: PARAMS['alfax_i']          = parameters['alfax_i']
        if 'alfay_i'          in parameters: PARAMS['alfay_i']          = parameters['alfay_i']
        if 'Ez'               in parameters: PARAMS['Ez_feld']          = parameters['Ez']
        if 'phi_sync'         in parameters: PARAMS['soll_phase']       = parameters['phi_sync']
        if 'sigmaz_i'         in parameters: PARAMS['sigmaz_i']         = parameters['sigmaz_i']
        if 'gap'              in parameters: PARAMS['spalt_laenge']     = parameters['gap']
        if 'cav_len'          in parameters: PARAMS['cavity_laenge']    = parameters['cav_len']
        if 'ql'               in parameters: PARAMS['ql']               = parameters['ql']
        if 'aperture'         in parameters: PARAMS['aperture']         = parameters['aperture']
        if 'windings'         in parameters: PARAMS['n_coil']           = parameters['windings']
        PARAMS['wellenlänge']    = PARAMS['lichtgeschwindigkeit']/PARAMS['frequenz']
        PARAMS['spalt_spannung'] = PARAMS['Ez_feld']*PARAMS['spalt_laenge']
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
            PARAMS['lattice_version'] = lattice_title
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
            PARAMS['nsuper'] = nsuper
            del segSubList[0]              #pull nsuper off
            # DEBUG('segSubList in expand_reduce()',segSubList)
            for i in range(nsuper):        #expand nsuper times
                for k in segSubList:
                    latticeList.append(k)
        return (latticeList,segments)
    ## factory body --------
    SUMMARY['input file'] = PARAMS['input_file'] = input_file

    with open(input_file,'r') as fileobject:
        in_data = yaml.load(fileobject)
    fileobject.close()

    read_flags(in_data)
    read_sections(in_data)
    read_parameters(in_data)
    # update PARAMS with overriding initials
    PARAMS.update(
        zellipse(PARAMS['sigmaz_i'],     ## calculate the long. emittance with def. parameters
                PARAMS['Ez_feld'],
                PARAMS['wellenlänge'],
        radians(PARAMS['soll_phase']),
                PARAMS['spalt_laenge'],
                PARAMS['sollteilchen']))
    PARAMS['emitz_i']  = PARAMS['emitz']   # here zellipse calculated initial values
    PARAMS['betaz_i']  = PARAMS['betaz']   # here zellipse calculated initial values
    PARAMS['gammaz_i'] = PARAMS['gammaz']  # here zellipse calculated initial values
    PARAMS['alfaz_i']  = PARAMS['alphaz']  # here zellipse calculated initial values
    # del unused key-values from zellipse
    del PARAMS['emitz']
    del PARAMS['betaz']
    del PARAMS['gammaz']
    del PARAMS['alphaz']
    # DEBUG('PARAMS after read_parameters()',PARAMS.__dict__)

    (latticeList,segments) = expand_reduce(in_data)
    # DEBUG('latticeList in factory()',latticeList)      # def of all segments in lattice
    # DEBUG('segments in factory()',segments)            # def of all segments

    PARAMS['sollteilchen'](tkin=PARAMS['injection_energy'])# (WICHTIG) set sollteilchen energy
    lattice = make_lattice(latticeList,segments)
    # DEBUG('lattice_generator >>\n',lattice.string())

    SUMMARY['aperture [m]']       = PARAMS['aperture']
    SUMMARY['lattice length [m]'] = PARAMS['lattice_length']  = lattice.length
    # DEBUG('SUMMARY in factory()',SUMMARY)
    return lattice    #end of factory(...)

def parse_yaml_and_fabric(input_file,factory=factory):   ## delegates to factory
    return factory(input_file)

## utilities
def test0():
    print('---------------------------------TEST0')
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
    print('---------------------------------TEST1')
    lattice = parse_yaml_and_fabric(input_file)
## main ----------
if __name__ == '__main__':
    test0()
    test1('fodo_with_10cav_per_RF(4).yml')

