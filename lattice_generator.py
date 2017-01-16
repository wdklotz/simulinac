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
from setup import CONF,SUMMARY,Particle,DEBUG
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
	# DEBUG('instanciate_element for item >>',item)
	key = item[0]
	att_list = item[1]
	attributes = unpack_list_of_dict(att_list)
# 	DEBUG('',key, attributes)
	if key == 'D':
		length   = attributes['length']
		label    = attributes['label']
		instance =  ELM.D(length=length,label=label,particle=Particle.soll)
	elif key == 'QF':
		length   = attributes['length']
		label    = attributes['label']
		dBdz     = attributes["B'"]
		kq       = dBdz/Particle.soll.brho
		instance = ELM.QF(k0=kq,length=length,label=label,particle=Particle.soll)
	elif key == 'QD':
		length   = attributes['length']
		label    = attributes['label']
		dBdz     = attributes["B'"]
		kq       = dBdz/Particle.soll.brho
		instance = ELM.QD(k0=kq,length=length,label=label,particle=Particle.soll)
	elif key == 'RFG':
		gap       = attributes['gap']
		label     = attributes['label']
		Ez        = attributes["Ez"]
		PhiSoll   = radians(attributes["PhiSync"])
		fRF       = attributes["fRF"]
		U0        = Ez * gap
		dWf       = CONF['dWf']
		instance  =  ELM.RFG(U0=U0,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,particle=Particle.soll,dWf=dWf)
	elif key == 'RFC':
		gap       = attributes['gap']
		length    = attributes['length']
		label     = attributes['label']
		Ez        = attributes["Ez"]
		PhiSoll   = radians(attributes["PhiSync"])
		fRF       = attributes["fRF"]
		U0        = Ez * gap
		dWf       = CONF['dWf']
		instance  =  ELM.RFC(U0=U0,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,length=length,particle=Particle.soll,dWf=dWf)
	else:
		raise RuntimeError('unknown element type: ',key)
	# DEBUG('{} instance created'.format(label))
	return (label,instance)

def factory(input_file):
#--------
	def make_lattice(ns,lat,seg,elm):
		lattice = Lattice()
		for h in range(ns):      #loop nsuper
			for i in lat:        #loop segments in lattice def
				seg_label = i
	# 			DEBUG('segment >>',seg_label)
				for j in seg:    #browse for segment def
					if j['label'] == seg_label:
						elm_list = j['elements']
						break
				for k in elm_list: #loop segment elements
					elm_label = k
					# DEBUG('\telement >>',elm_label)
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
# --------
# --------
	def read_flags(in_data):
	#returns ==> {...}
		flags_list = in_data['flags']
		flags      = unpack_list_of_dict(flags_list)
	# 	DEBUG('flags=\t',flags)
		CONF['dWf'] = SUMMARY['acc. ON']                   = flags['accON']
		CONF['periodic'] = SUMMARY['ring lattice']         = flags['periodic']
		CONF['verbose']                                    = flags['verbose']
		return flags
# --------
	def read_parameters(in_data):
	#returns ==> {...}
		parameter_list = in_data['parameters']
		parameters     = unpack_list_of_dict(parameter_list)
		# DEBUG('parameters=\t',parameters)
		CONF['frequenz']         = parameters['frequency']
		CONF['quad_gradient']    = None if not 'B_grad' in parameters else parameters['B_grad']
		CONF['quadf_gradient']   = CONF['quad_gradient'] if not 'B_grad_f' in parameters else parameters['B_grad_f']
		CONF['quadd_gradient']   = CONF['quad_gradient'] if not 'B_grad_d' in parameters else parameters['B_grad_d']
		CONF['injection_energy'] = parameters['Tkin']
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
# --------
	def expand_reduce(in_data):
	#--------
		def read_elements(in_data):
			elements_list = in_data['elements']
			elements_dict = unpack_list_of_dict(elements_list)
			return elements_dict
	#--------
		def read_segments(in_data):
			segments_list = in_data['segments']
			segments_dict = unpack_list_of_dict(segments_list)
			return segments_dict
	#--------
		def read_lattice(in_data):
			lattice_segment_list= in_data['lattice']
			lattice_title = lattice_segment_list[0]['label']
			del lattice_segment_list[0]         #pull label off
			CONF['lattice_version'] = lattice_title
			return lattice_segment_list[0]
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
	# 				DEBUG('e_list_item ==> ',e_list_item)
					e_label = e_list_item[0]['label']
					list_of_segment_items.append(e_label)
				segment['elements']=list_of_segment_items
				list_of_segments.append(segment)
			return list_of_segments
	#--------
		lattice_segments = read_lattice(in_data)
		nsuper = lattice_segments[0]       #nboff super cells
		del lattice_segments[0]            #pull nsuper off
		seg_defs = read_segments(in_data)
		segments = reduce_seg_def(seg_defs)
		elem_defs = read_elements(in_data)
		elements = reduce_elm_def(elem_defs)
		return (nsuper,lattice_segments,segments,elements)
	#--------
	SUMMARY['input file'] = CONF['input_file'] = input_file

	with open(input_file,'r') as fileobject:
		in_data = yaml.load(fileobject)
	fileobject.close()

	read_flags(in_data)
	read_parameters(in_data)
	(nsuper,lattice_in,segments_in, elements_in) = expand_reduce(in_data)
# 	DEBUG('\nnsuper ==>',nsuper)   #nboff super cells
# 	DEBUG('\nlattice_segments ==>',lattice_in)  #def of all segments in lattice
# 	DEBUG('\nsegments ==>',segments_in)  #def of all segments
# 	DEBUG('\nelements ==>',elements_in)  #def of all elements
	lattice = make_lattice(nsuper,lattice_in,segments_in,elements_in)
# 	DEBUG('lattice >>',end='\n')
	lattice.string()
	SUMMARY['lattice length [m]'] = CONF['lattice_length']  = lattice.length
	return lattice    #end of factory(...)

def parse_yaml_and_fabric(input_file,factory=factory):   ## delegates to factory
	return factory(input_file)

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
#--------
if __name__ == '__main__':
#     test0()
	test1('fodo_with_10cav_per_RF(2).yml')

