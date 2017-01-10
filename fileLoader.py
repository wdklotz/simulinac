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
from setup import CONF,SUMMARY,Beam
from setup import objprnt,Wakzeptanz
import elements as ELM
from lattice import Lattice
from math import radians,sqrt,pi,degrees
import yaml

def unpack_list_of_dict(alist):
	new = {}
	for item in alist:
		for key in list(item.keys()):
			new[key] = item[key]
	return new

def instanciate_element(item):
# 	print('item >>',item)
	key = item[0]
	att_list = item[1]
	attributes = unpack_list_of_dict(att_list)
# 	print(key, attributes)
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
		segment_id = item[0]
# 		print(segment_id)
		element_list = item[1]
# 		print(element_list)
		segment_label = element_list[0]['label']        ## pull {'label:xxx'} off
# 		print(segment_label)
		del element_list[0]
# 		print(element_list)
		segment = Lattice()
		for element in element_list:
			attributes = unpack_list_of_dict(element)
# 			print(attributes)
			instance = instances_dict[attributes['label']]
			segment.add_element(instance)
			segment_instance_dict[segment_label] = segment
	return segment_instance_dict

def read_flags(in_data):
#returns ==> {...}
	flags_list = in_data['flags']
	flags      = unpack_list_of_dict(flags_list)
# 	print('\nflags=\t',flags)
	CONF['dWf'] = SUMMARY['acc. ON']                   = flags['accON']
	CONF['periodic'] = SUMMARY['ring lattice']         = flags['periodic']
	CONF['verbose']                                    = flags['verbose']
	return flags

def read_parameters(in_data):
#returns ==> {...}
	parameter_list = in_data['parameters']
	parameters     = unpack_list_of_dict(parameter_list)
	# print('parameters=\t',parameters)
	CONF['frequenz']         = parameters['frequency']
	CONF['quad_gradient']    = None if not 'B_grad' in parameters else parameters['B_grad']
	CONF['quadf_gradient']   = CONF['quad_gradient'] if not 'B_grad_f' in parameters else parameters['B_grad_f']
	CONF['quadd_gradient']   = CONF['quad_gradient'] if not 'B_grad_d' in parameters else parameters['B_grad_d']
	CONF['injection_energy'] = parameters['TK_i']
	CONF['emitx_i']          = parameters['emitx_i']
	CONF['emity_i']          = parameters['emity_i']
	CONF['betax_i']          = parameters['betax_i']
	CONF['betay_i']          = parameters['betay_i']
	CONF["alfax_i"]          = parameters['alfax_i']
	CONF["alfay_i"]          = parameters['alfay_i']
	CONF['dP/P']             = parameters['dP/P'] * 1.e-2
	CONF['Ez_feld']          = parameters['Ez']
	CONF['soll_phase']       = parameters['phi_sync']
	CONF['dZ']               = parameters['dZ']
	CONF['spalt_laenge']     = parameters['gap']
	CONF['cavity_laenge']    = parameters['cav_len']
	CONF['ql']               = parameters['ql']
	CONF['wellenlänge']      = CONF['lichtgeschwindigkeit']/CONF['frequenz']
	CONF['spalt_spannung']   = CONF['Ez_feld']*CONF['spalt_laenge']
	CONF['n_coil']           = 1 if not 'windings' in parameters else parameters['windings']
	return parameters

def read_elements(in_data):
	elements_list = in_data['elements']
	elements_dict = unpack_list_of_dict(elements_list)
	return elements_dict

def read_segments(in_data):
	segments_list = in_data['segments']
	segments_dict = unpack_list_of_dict(segments_list)
	return segments_dict

def read_lattice(in_data):
	lattice_segment_list= in_data['lattice']
	lattice_title = lattice_segment_list[0]['label']
	del lattice_segment_list[0]         #pull label off
	CONF['lattice_version'] = lattice_title
	return lattice_segment_list[0]

def expand_reduce(in_data):
#--------
	def merge_list_of_dicts(lstofdicts):
	#returns ==> {...}
		res={}
		for d in lstofdicts:
			for k,v in d.items():
				res[k]=v
		return res
#--------
	def reduce_elm_def(dict):
	#returns ==> [{k:v,...,k:v},{k:v,...,k:v},...,{k:v,...,k:v}]
		res=[]
		for k,v in dict.items():
			type = k
			p_list = v
			p_dict = merge_list_of_dicts(p_list)
			p_dict['type'] = type
			res.append(p_dict)
		return res
#--------
	def reduce_seg_def(dict):
	#returns ==> [{"label":str,"elements":[.....]},...,{"label":str,"elements":[.....]}]
		list_of_segments=[]
		for key,l in dict.items():
			segment={}
			outer_list = l
			seg_label = outer_list[0]['label']
			segment['label'] = seg_label
			e_list = outer_list[1:]
			list_of_segment_items=[]
			for e_list_item in e_list:
# 				print('e_list_item ==> ',e_list_item)
				e_label = e_list_item[0]['label']
				list_of_segment_items.append(e_label)
			segment['elements']=list_of_segment_items
			list_of_segments.append(segment)
		return list_of_segments
#--------
	lattice_segments = read_lattice(in_data)
	nsuper = lattice_segments[0]       #nboff super cells
	del lattice_segments[0]         #pull nsuper off
	seg_defs = read_segments(in_data)
	segments = reduce_seg_def(seg_defs)
	elem_defs = read_elements(in_data)
	elements = reduce_elm_def(elem_defs)
	return (nsuper,lattice_segments,segments,elements)

def collect_summaries():
	SUMMARY['frequency [Hz]'] = CONF['frequenz']
	SUMMARY['Quad gradient [T/m]'] = CONF['quad_gradient']
	SUMMARY['QF gradient [T/m]'] = CONF['quadf_gradient']
	SUMMARY['QD gradient [T/m]'] = CONF['quadd_gradient']
	SUMMARY['Quad pole length [m]'] = CONF['ql']
	SUMMARY['injection energy [MeV]'] = CONF['injection_energy']
	SUMMARY['emitx_i [rad*m]'] = CONF['emitx_i']
	SUMMARY['emity_i [rad*m]'] = CONF['emity_i']
	SUMMARY['sigx_i* [mm]'] = 1000.*sqrt(CONF['betax_i']*CONF['emitx_i'])  # enveloppe @ entrance
	SUMMARY['sigy_i* [mm]'] = 1000.*sqrt(CONF['betay_i']*CONF['emity_i'])
	SUMMARY['<impulse dP/P> [%]'] = CONF['dP/P'] * 1.e+2
	SUMMARY['sync. phase [deg]'] = CONF['soll_phase']
	SUMMARY['dZ [m]'] = CONF['dZ']
	SUMMARY['cavity gap length [m]'] = CONF['spalt_laenge']
	SUMMARY['cavity length [m]'] = CONF['cavity_laenge']
	SUMMARY['wavelength [m]'] = CONF['wellenlänge']
	SUMMARY['cavity gap voltage* [MV]'] = CONF['spalt_spannung']
	SUMMARY['acc. field Ez [MV/m]'] = CONF['Ez_feld']
	SUMMARY['lattice version'] = CONF['lattice_version']
	#...........*...........*...........*...........*...........*...........*
	SUMMARY['QF pole strength* [T]'] = CONF['quadf_gradient'] * CONF['ql']
	SUMMARY['QF current* [A/winding]'] = (CONF['quadf_gradient'] * (CONF['ql']*1000.)**2 )/2.52/CONF['n_coil']
	SUMMARY['QF power estimate* [W]'] = 0.0115 *SUMMARY['QF current* [A/winding]']**2  # R=0.0115 Ohms
	SUMMARY['QF coil [windings]'] = CONF['n_coil']
	SUMMARY['<energy dW/W> max* [%]'] = wakzp = Wakzeptanz(    # energy acceptance in %
		CONF['Ez_feld'],
		Beam.soll.TrTf(CONF['spalt_laenge'],CONF['frequenz']),
		CONF['soll_phase'],
		CONF['wellenlänge'],
		Beam.soll)*1.e+2
	SUMMARY['<impulse dP/P> max* [%]'] = 1./(1.+1./Beam.soll.gamma)*wakzp  # impule acceptanc in %
	SUMMARY['dphi* [deg]'] =degrees( 2*pi*CONF['frequenz']/CONF['lichtgeschwindigkeit']/Beam.soll.beta*CONF['dZ'])

# 	dictprnt(SUMMARY,text='summary')
	return

def parse_yaml_and_fabric(factory,input_file):          ## delegates
	return factory(input_file)

def factory1(input_file):          ## lattice factory v1
#--------
	def make_lattice(lattice_segment_list,segment_instance_dict):
		lattice = Lattice()
		nsuper = lattice_segment_list[0]           ## pull nsuper off...
		del lattice_segment_list[0]
		seg_counter = 0
		for anz in range(nsuper):
			for segment_label in lattice_segment_list:
				lattice_part = segment_instance_dict[segment_label]
				lattice.append(lattice_part)
				seg_counter += 1
		SUMMARY['nboff segments*']= seg_counter
		lattice.out()
		return lattice
#--------
	SUMMARY['input file'] = CONF['input_file'] = input_file
	with open(input_file,'r') as fileobject:
		in_data = yaml.load(fileobject)
#--------
	flags        = read_flags(in_data)
	parameters   = read_parameters(in_data)
#--------
	elements_dict = read_elements(in_data)
#     print('\nelements=\t',elements_dict)
	instances_dict = {}
	for item in elements_dict.items():
		(label,instance) = instanciate_element(item)
		instances_dict[label]= instance
#     print('\ninstances=\t',instances_dict.keys())
#--------
	segments_dict = read_segments(in_data)
#     print('\nsegments=\t',segments_dict)
	segment_instance_dict= make_segments(segments_dict,instances_dict)
#     print('segment_instances=\t',segment_instance_dict)
#--------
    # proton: the default synchronous reference particle  (class member!)
# 	Beam.soll             = Proton(CONF['injection_energy'])
    # objprnt(Beam.soll,text='injected beam')
#--------
	lattice_segment_list = read_lattice(in_data)
	lattice = make_lattice(lattice_segment_list,segment_instance_dict)
	collect_summaries()
	return lattice

def factory2(input_file):
#--------
	def make_lattice(ns,lat,seg,elm):
		lattice = Lattice()
		for h in range(ns):      #loop nsuper
			for i in lat:        #loop segments in lattice def
				seg_label = i
	# 			print('segment >>',seg_label)
				for j in seg:    #browse for segment def
					if j['label'] == seg_label:
						elm_list = j['elements']
						break
				for k in elm_list: #loop segment elements
					elm_label = k
	# 				print('\telement >>',elm_label)
					for l in elm: #browse for element in seg
						if l['label'] == elm_label:
							attributes=[]
							for m,n in l.items():  #build item description
								attributes.append({m:n})
							item = (l['type'],attributes)
							label,instance = instanciate_element(item)  #instanciate
							lattice.add_element(instance)  #add element instance to lattice
							break
		return lattice   #the complete lattice
#--------
	SUMMARY['input file'] = CONF['input_file'] = input_file

	with open(input_file,'r') as fileobject:
		in_data = yaml.load(fileobject)
	fileobject.close()

	read_flags(in_data)
	read_parameters(in_data)
	(nsuper,lattice_in,segments_in, elements_in) = expand_reduce(in_data)
# 	print('nsuper ==>',nsuper)   #nboff super cells
# 	print('\nlattice_segments ==>',lattice_in)  #def of all segments in lattice
# 	print('\nsegments ==>',segments_in)  #def of all segments
# 	print('\nelements ==>',elements_in)  #def of all elements
	lattice = make_lattice(nsuper,lattice_in,segments_in,elements_in)
# 	print('lattice >>',end='\n')
	lattice.out()
	collect_summaries()
	SUMMARY['lattice length [m]'] = CONF['lattice_length']  = lattice.length
	return lattice

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
	lattice = parse_yaml_and_fabric(factory2,input_file)

def test2(input_file):
	print('\nTest2')
	parse_yaml_and_fabric(factory2,input_file)
#--------
if __name__ == '__main__':
#     test0()
# 	test1('fodo_with_10cav_per_RF(2).yml')
	test2('fodo_with_10cav_per_RF(2).yml')

