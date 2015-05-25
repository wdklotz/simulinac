#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import CONF,Beam,Proton,objprnt
import elements as ELM
from lattice import Lattice
import yaml
from math import radians

def unpack_list_of_dict(alist):
    new = {}
    for item in alist:
        for key in list(item.keys()):
            new[key] = item[key]
    return new
def instanciate_element(item):
    key = item[0]
    att_list = item[1]
    attributes = unpack_list_of_dict(att_list)
    # print(key, attributes)
    if key == 'D':
        length   = attributes['length']
        label    = attributes['label']
        instance =  ELM.D(length=length,label=label,beam=Beam.soll)
        return (label,instance)
    if key == 'QF':
        length   = attributes['length']
        label    = attributes['label']
        dBdz     = attributes["B'"]
        kq       = dBdz/Beam.soll.brho
        instance = ELM.QF(k0=kq,length=length,label=label,beam=Beam.soll)
        return (label,instance)
    if key == 'QD':
        length   = attributes['length']
        label    = attributes['label']
        dBdz     = attributes["B'"]
        kq       = dBdz/Beam.soll.brho
        instance = ELM.QD(k0=kq,length=length,label=label,beam=Beam.soll)
        return (label,instance)
    if key == 'RFG':
        gap       = attributes['gap']
        label     = attributes['label']
        Ez        = attributes["Ez"]
        PhiSoll   = radians(attributes["PhiSync"])
        fRF       = attributes["fRF"]
        U0        = Ez * gap
        dWf       = CONF['dWf']
        instance  =  ELM.RFG(U0=U0,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,beam=Beam.soll,dWf=dWf)
        return (label,instance)
    if key == 'RFC':
        gap       = attributes['gap']
        length    = attributes['length']
        label     = attributes['label']
        Ez        = attributes["Ez"]
        PhiSoll   = radians(attributes["PhiSync"])
        fRF       = attributes["fRF"]
        U0        = Ez * gap
        dWf       = CONF['dWf']
        instance  =  ELM.RFC(U0=U0,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,length=length,beam=Beam.soll,dWf=dWf)
        return (label,instance)
    else:
        raise RuntimeError('unknown element type: ',key)
def make_segments(segments_dict,instances_dict):
    segment_instance_dict = {}
    for item in segments_dict.items():
        # segment_id = item[0]
        # print(segment_id)
        element_list = item[1]
        # print(element_list)
        segment_label = element_list[0]['label']
        # print(segment_label)
        del element_list[0]
        # print(element_list)
        segment = Lattice()
        for element in element_list:
            attributes = unpack_list_of_dict(element)
            # print(attributes)
            instance = instances_dict[attributes['label']]
            segment.add_element(instance)
        segment_instance_dict[segment_label] = segment
    return segment_instance_dict
def make_lattice(lattice_segment_list,segment_instance_dict):
    lattice = Lattice()
    for segment_label in lattice_segment_list:
        lattice_part= segment_instance_dict[segment_label]
        lattice.append(lattice_part)
    return lattice
def test0():
    wfl= []
    fileobject=open('template.yml','r')
    wfl= yaml.load(fileobject)
    print(yaml.dump(wfl,default_flow_style=True))
    for i,v in iter(wfl.items()):
        print(i,' =\t',v)
    seg = wfl['segment']
    print(seg)
    print('=== segment ===')
    for i in seg:
        print(i)
    lattice = wfl['lattice']
    print('=== lattice ===')
    for l in lattice:
        for i in l:
            print(i)
def read_yaml_and_parse(filepath):
    fileobject = open(filepath,'r')
    in_data    = yaml.load(fileobject)
#...........*...........*...........*...........*...........*...........*...........*
    flags_list = in_data['flags']
    flags      = unpack_list_of_dict(flags_list)
    # print('\nflags=\t',flags)
    CONF['dWf']              = flags['accON']
    CONF['periodic']         = flags['periodic']
#...........*...........*...........*...........*...........*...........*...........*
    parameter_list = in_data['parameters']
    parameters     = unpack_list_of_dict(parameter_list)
    # print('parameters=\t',parameters)
    CONF['frequenz']         = parameters['frequency']
    CONF['quad_gradient']    = parameters['B_grad']
    CONF['injection_energy'] = parameters['TK_i']
    CONF['emitx_i']          = parameters['emitx_i']
    CONF['emity_i']          = parameters['emity_i']
    CONF['sigx_i']           = parameters['sigx_i']
    CONF['sigy_i']           = parameters['sigy_i']
    CONF['dP/P']             = parameters['dP/P'] * 1.e-2
    CONF['Ez_feld']          = parameters['Ez']
    CONF['soll_phase']       = parameters['phi_sync']
    CONF['dZ']               = parameters['dZ']
    CONF['spalt_laenge']     = parameters['gap']
    CONF['cavity_laenge']    = parameters['cav_len']
    
    CONF['emity_i']       = CONF['emitx_i']
    CONF['sigy_i']        = CONF['sigx_i']
    CONF['wellenlänge']   = CONF['lichtgeschwindigkeit']/CONF['frequenz']
    CONF['spalt_spannung']= CONF['Ez_feld']*CONF['spalt_laenge']
#...........*...........*...........*...........*...........*...........*...........*
    # proton: the default synchronous reference particle  (class member!)
    Beam.soll             = Proton(CONF['injection_energy'])
    # objprnt(Beam.soll,text='injected beam')
#...........*...........*...........*...........*...........*...........*...........*
    elements_list = in_data['elements']
    elements_dict = unpack_list_of_dict(elements_list)
    # print('\nelements=\t',elements_dict)
    instances_dict = {}
    for item in elements_dict.items():
        (label,instance) = instanciate_element(item)
        instances_dict[label]= instance
    # print(instances_dict.keys())
#...........*...........*...........*...........*...........*...........*...........*
    segments_list = in_data['segments']
    segments_dict = unpack_list_of_dict(segments_list)
    # print('\nsegments=\t',segments_dict)
    segment_instance_dict= make_segments(segments_dict,instances_dict)
    # print(segment_instance_dict)
#...........*...........*...........*...........*...........*...........*...........*
    lattice_segment_list= in_data['lattice']
    # print('\nlattice=\t',lattice_segment_list)
    lattice_title = lattice_segment_list[0]['label']
    del lattice_segment_list[0]
    # print('\nlattice=\t',lattice_segment_list)
    lattice = make_lattice(lattice_segment_list,segment_instance_dict)
    lattice.trim_energy()
    print(lattice_title)
    # lattice.out()
    return lattice
#...........*...........*...........*...........*...........*...........*...........*
if __name__ == '__main__':
#     test0()
    lattice = read_yaml_and_parse('template.yml')
    lattice.out()

