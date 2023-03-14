#!/Users/klotz/SIMULINAC_env/bin/python
# -*- coding: utf-8 -*-
__version__='v10.23.2'
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
from sre_parse import FLAGS
import sys,os
from math import radians
import yaml
import warnings
import unittest

import setutil             as UTIL
import elements            as ELM
import OXAL                as OXA
import TTFG                as TTF
import DYNG                as DYN
import PsMarkerAgent       as PSMKR
import PoincareMarkerAgent as PCMKR
from lattice import Lattice
from Ez0 import SFdata
from lattice_parser2 import parse as doInputParser

def make_counter():
    count = 0
    def inner():
        nonlocal count
        count += 1
        return count
    return inner
counter1 = make_counter()
counter2 = make_counter()

def marker_is_compatible_with(prog,ID):
    ret = True
    head, tail = os.path.split(sys.argv[0])
    this_prog = tail
    UTIL.DEBUG_OFF((prog,this_prog))
    if prog != this_prog:
        ret = False
    return ret
def get_mandatory(attributes,key,item):
    try:
        res = attributes[key]
    except KeyError:
        warnings.showwarning(
                'InputError: Mandatory attribute "{}" missing for element "{}" - STOP'.format(key,item),
                UserWarning,
                'lattice_generator.py',
                'get_mandatory()',
                )
        sys.exit(1)
    return res
def instanciate_element(item):
    """ item: {ID:{attributes}} for each node """
    instance = None     # will be defined below and returned
    for ID,attributes in item.items():
        UTIL.DEBUG_OFF(F"ID={ID} attributes={attributes}")
        ELEMENT = UTIL.ELEMENTS[ID]          # the item in the ELEMENT list
        type = attributes.get('type')
        if type   == 'D':
            length         = get_mandatory(attributes,'length',ID)
            aperture       = attributes.get('aperture')
            instance       =  ELM.D(ID,length=length,aperture=aperture)
            ELEMENT['sec'] = attributes.get('sec','?')
        elif type == 'DKD':
            length         = get_mandatory(attributes,'length',ID)
            aperture       = attributes.get('aperture')
            instance       =  ELM.DKD(ID,length=length,aperture=aperture)
            ELEMENT['sec'] = attributes.get('sec','?')
        elif type == 'QF':
            length           = get_mandatory(attributes,'length',ID)
            dBdz             = get_mandatory(attributes,"B'",ID)
            aperture         = get_mandatory(attributes,'aperture',ID)
            instance         = ELM.QF(ID,dBdz,length=length,aperture=aperture)
            ELEMENT['Bpole'] = dBdz*aperture      # Bpole
            ELEMENT['sec']   = attributes.get('sec','?')
        elif type == 'QD':
            length           = get_mandatory(attributes,'length',ID)
            dBdz             = get_mandatory(attributes,"B'",ID)
            aperture         = get_mandatory(attributes,'aperture',ID)
            instance         = ELM.QD(ID,dBdz,length=length,aperture=aperture)
            ELEMENT['Bpole'] = dBdz*aperture      # Bpole
            ELEMENT['sec']   = attributes.get('sec','?')
        elif type == 'RFG':
            phiSoll   = radians(get_mandatory(attributes,"PhiSync",ID))
            freq      = float(get_mandatory(attributes,"freq",ID))
            gap       = get_mandatory(attributes,'gap',ID)
            aperture  = get_mandatory(attributes,'aperture',ID)
            dWf       = UTIL.FLAGS['dWf']
            EzPeak    = get_mandatory(attributes,"EzPeak",ID)
            mapping   = UTIL.FLAGS.get('mapping')    # global mapping FLAG overrides individual mapping
            if mapping == None:
                mapping = attributes.get('mapping','t3d')
            UTIL.ELEMENTS[ID]['mapping'] = mapping   # maybe overriden by global mapping
            if mapping == 'ttf' or mapping == 'dyn' or mapping == 'oxal': # SF-data
                UTIL.FLAGS['non_linear_mapping'] = True
                fieldtab = get_mandatory(attributes,"SFdata",ID)
                gap_cm = gap*100     # Watch out!
                sfdata = SFdata.field_data(fieldtab,EzPeak=EzPeak,gap=gap_cm)
            else:
                ELEMENT['EzAvg']  = EzAvg = EzPeak
                ELEMENT['SFdata'] = None
            if mapping   == 'oxal':
                EzAvg = ELEMENT['EzAvg'] = sfdata.EzAvg
                instance = OXA.OXAL_G(ID,EzAvg,phiSoll,gap,freq,SFdata=sfdata,particle=UTIL.Proton(UTIL.PARAMS['injection_energy']),position=(0.,0.,0.),aperture=aperture,dWf=dWf,fieldtab=fieldtab)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            elif mapping == 'ttf':
                EzAvg = ELEMENT['EzAvg'] = sfdata.EzAvg
                instance = TTF.TTF_G(ID,EzAvg,phiSoll,gap,freq,SFdata=sfdata,particle=UTIL.Proton(UTIL.PARAMS['injection_energy']),position=(0.,0.,0.),aperture=aperture,dWf=dWf,fieldtab=fieldtab)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            elif mapping == 'dyn':
                if counter1() < 5:
                    print(UTIL.colors.RED+"WARN: dyn mapping is broken"+UTIL.colors.ENDC)
                EzAvg = ELEMENT['EzAvg'] = sfdata.EzAvg
                instance = DYN.DYN_G(ID,EzAvg,phiSoll,gap,freq,SFdata=sfdata,particle=UTIL.Proton(UTIL.PARAMS['injection_energy']),position=(0.,0.,0.),aperture=aperture,dWf=dWf)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            else:
                instance = ELM.RFG(ID,EzAvg,phiSoll,gap,freq,SFdata=None,particle=UTIL.Proton(UTIL.PARAMS['injection_energy']),position=(0.,0.,0.),aperture=aperture,dWf=dWf,mapping=mapping)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            """         elif type == 'RFC':
            if counter2() < 5:
                print(UTIL.colors.RED+"WARN: RFC node is broken"+UTIL.colors.ENDC)
            else:
                sys.exit(1)
            phiSoll   = radians(get_mandatory(attributes,"PhiSync",ID))
            freq      = float(get_mandatory(attributes,"freq",ID))
            gap       = get_mandatory(attributes,'gap',ID)
            length    = get_mandatory(attributes,'length',ID)
            aperture  = get_mandatory(attributes,'aperture',ID)
            dWf       = UTIL.FLAGS['dWf']
            EzPeak    = get_mandatory(attributes,"EzPeak",ID)
            mapping   = attributes.get('mapping','t3d')
            if mapping == 'ttf' or mapping == 'dyn' or mapping == 'oxal': # SF-data
                UTIL.FLAGS['non_linear_mapping'] = True
                fieldtab = get_mandatory(attributes,"SFdata",ID)
                gap_cm = gap*100     # Watch out!
                sfdata = SFdata.field_data(fieldtab,EzPeak=EzPeak,gap=gap_cm)
                # if fieldtab not in UTIL.PARAMS:
                #     gap_cm = gap*100     # Watch out!
                #     UTIL.PARAMS[fieldtab] = SFdata(fieldtab,EzPeak=EzPeak,gap=gap_cm)
            else:
                ELEMENT['EzAvg']  = EzAvg = EzPeak
                ELEMENT['SFdata'] = None
            if mapping == 'oxal':
                EzAvg = ELEMENT['EzAvg'] = UTIL.PARAMS[fieldtab].EzAvg
                # instance = OXA.OXAL_C(ID,EzAvg,phiSoll,gap,length,freq,SFdata=UTIL.PARAMS[fieldtab],aperture=aperture,dWf=dWf,fieldtab=fieldtab)
                instance = OXA.OXAL_C(ID,EzAvg,phiSoll,gap,length,freq,SFdata=sfdata,aperture=aperture,dWf=dWf,fieldtab=fieldtab)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            elif mapping == 'ttf':
                EzAvg = ELEMENT['EzAvg'] = UTIL.PARAMS[fieldtab].EzAvg
                # instance = TTF.TTF_C(ID,EzAvg,phiSoll,gap,length,freq,SFdata=UTIL.PARAMS[fieldtab],aperture=aperture,dWf=dWf,fieldtab=fieldtab)
                instance = TTF.TTF_C(ID,EzAvg,phiSoll,gap,length,freq,SFdata=sfdata,aperture=aperture,dWf=dWf,fieldtab=fieldtab)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            elif mapping == 'dyn':
                EzAvg = ELEMENT['EzAvg'] = UTIL.PARAMS[fieldtab].EzAvg
                # instance = DYN.DYN_C(ID,EzAvg,phiSoll,gap,length,freq,SFdata=UTIL.PARAMS[fieldtab],aperture=aperture,dWf=dWf,fieldtab=fieldtab)
                instance = DYN.DYN_C(ID,EzAvg,phiSoll,gap,length,freq,SFdata=sfdata,aperture=aperture,dWf=dWf,fieldtab=fieldtab)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            else:
                instance = ELM.RFC(ID,EzAvg,phiSoll,gap,freq,length,mapping=mapping,aperture=aperture,dWf=dWf,fieldtab=None)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            """        
        elif type == 'GAP':
            gap       = get_mandatory(attributes,'gap',ID)
            EzPeak    = get_mandatory(attributes,"EzPeak",ID)
            phiSoll   = radians(get_mandatory(attributes,"PhiSync",ID))
            freq      = float(get_mandatory(attributes,"freq",ID))
            dWf       = UTIL.FLAGS['dWf']
            aperture  = get_mandatory(attributes,'aperture',ID)
            EzAvg     = EzPeak
            instance  =  ELM.GAP(ID,EzAvg,phiSoll,gap,freq,aperture=aperture,dWf=dWf)
            ELEMENT['EzPeak'] = EzPeak
            ELEMENT['sec']    = attributes.get('sec','?')
            ELEMENT['EzAvg']  = EzAvg
        elif type == 'MRK':
            active  = attributes.get('active',UTIL.FLAGS['mactive'])
            viseo   = attributes.get('viseo',3)
            ELEMENT = UTIL.ELEMENTS[ID]
            action  = get_mandatory(attributes,'action',ID)
            if action == 'pspace':
                # A marker for simu.py ?
                if not marker_is_compatible_with('simu.py',ID):
                    active = False
                    UTIL.DEBUG_OFF(UTIL.colors.RED+f'WARN: Marker {ID} incompatible with simu.py. Will be skipped'+UTIL.colors.ENDC)
                instance = PSMKR.PsMarkerAgent(ID,active,viseo)
                sec = attributes.get('sec','?') 
                ELEMENT['sec']   = sec
                UTIL.DEBUG_OFF(ELEMENT)
                UTIL.DEBUG_OFF(instance.toString())

            elif action == 'pcrcut':
                # A marker for tracker.py ?
                if not marker_is_compatible_with('tracker.py',ID):
                    active = False
                    UTIL.DEBUG_OFF(UTIL.colors.RED+f'WARN: Marker {ID} incompatible with tracker.py. Will be skipped'+UTIL.colors.ENDC)
                sec        = attributes.get('sec','?') 
                prefix     = attributes.get('prefix','frames')
                abscissa   = attributes.get('abscissa','z')
                ordinate   = attributes.get('ordinate','zp')
                instance   = PCMKR.PoincareMarkerAgent(ID,active,viseo,prefix,abscissa,ordinate)
                ELEMENT['sec']      = sec
                ELEMENT['prefix']   = prefix
                ELEMENT['abscissa'] = abscissa
                ELEMENT['ordinate'] = ordinate
                UTIL.DEBUG_OFF(ELEMENT)
                UTIL.DEBUG_OFF(instance.__dict__)
            else:
                warnings.showwarning(
                        'InputError: Unknown marker ACTION encountered: "{}" - STOP'.format(action),
                        UserWarning,
                        'lattice_generator.py',
                        'instanciate_element()',
                        )
                sys.exit(1)
        else:
            warnings.showwarning(
                'InputError: Unknown element TYPE encountered: "{}" - STOP'.format(type),
                UserWarning,
                'lattice_generator.py',
                'instanciate_element()',
                )
            sys.exit(1)
    return instance
def factory(input_file,stop=None):
    """ factory reads FLAGS, PARAMETERS and creates a lattice from yaml input """
    def proces_flags(flags):
        """ external FLAGs """        
        res = dict(
            accON    = flags.get('accON',True),               # acceleration ON
            periodic = flags.get('periodic',False),           # periodic lattice? default
            egf      = flags.get('egf',False),                # emittance grow flag default
            sigma    = flags.get('sigma',True),               # beam sizes by sigma-tracking
            KVout    = flags.get('KVout',False),              # print a dictionary of Key-Value pairs, no display
            verbose  = flags.get('verbose',0),                # print flag default = 0
            useaper  = flags.get('useaper',False),            # use aperture check for quads and rf-gaps
            bucket   = flags.get('bucket',False),             # plot bucket
            csTrak   = flags.get('csTrak',True),              # plot CS trajectories
            mactive  = flags.get('mactive',False),            # call marker actions
            envelope = flags.get('envelope',False),           # plot transverse envelopes
            mapping  = flags.get('mapping'),                  # global mapping overrides individula mapping (default to None)
        )
        """ internal FLAGs """        
        res['non_linear_mapping'] = False
        return res
    def proces_parameters(parameters):
        """ external parameters and their default values """
        res = dict()
        # kinetic energy @ injection
        W_in  = parameters.get("Win",UTIL.PARAMS['injection_energy']) # default PARAMS['injection_energy']
        tk_in = parameters.get('Tkin',W_in) # alias Tkin=Win
        res['injection_energy'] = tk_in 
        # longitudinal in {Dphi,w}-space
        DW2W_in = parameters.get('DW2W',0.01) # default 1%
        DT2T_in = parameters.get('DT2T',DW2W_in) # alias DT2T=DW2W
        res['DT2T']             = DT2T_in
        res['Dphi0']            = radians(parameters.get('DPHI0',10.)) # default [rad]
        # transverse beam parameters
        res['emitx_i']          = parameters.get('emitx_i',10E-6) # [m*rad]
        res['betax_i']          = parameters.get('betax_i',1.) # [m]
        res['alfax_i']          = parameters.get('alfax_i',0.)
        res['emity_i']          = parameters.get('emity_i',10E-6)
        res['betay_i']          = parameters.get('betay_i',1.)
        res['alfay_i']          = parameters.get('alfay_i',0.)
        # transverse Twiss @ entrance
        res['twiss_x_i'] = UTIL.Twiss(res['betax_i'], res['alfax_i'],res['emitx_i'])
        res['twiss_y_i'] = UTIL.Twiss(res['betay_i'], res['alfay_i'],res['emity_i'])
        # supplemental global parameters
        res['nbsigma']          = parameters.get('nbsigma',2)
        res['lattice_version']  = parameters.get('lattvers','not given')
        res['thins']            = parameters.get('thins',1)
        res['input_file']       = None
        # longitudinal emittance @ entrance
        Dphi0   = res['Dphi0']
        DT2T    = res['DT2T']
        T       = UTIL.PARAMS['injection_energy']
        E0      = UTIL.PARAMS['proton_mass']
        w0      = T/E0*DT2T # Wrangler's definition of w (pp.176)
        emit_w  = Dphi0 * w0 # emittance in {Dphi,w}-space
        res['w0']       = w0
        res['emitw_i']  = Dphi0 * w0
        res['alfaw_i']  = 0.   # always
        res['betaw_i']  = emit_w/w0**2
        # longitudinal TWiss @ entrance
        res['twiss_w_i'] = UTIL.Twiss(res['betaw_i'], res['alfaw_i'],res['emitw_i'])
        # set parameters that will be caculated later (f.i. waccept) to None
        res['emitz_i']   = None
        res['alfaz_i']   = None
        res['betaz_i']   = None
        return res
    def proces_elements(elements):
        """fills global ELEMENTS"""
        return elements
    def make_lattice(elementIDs):
        UTIL.DEBUG_OFF(elementIDs)
        lattice = Lattice(descriptor=UTIL.PARAMS.get('descriptor'))
        instances = []
        for elementID in elementIDs:
            UTIL.DEBUG_OFF(elementID)
            ELEMENT = UTIL.ELEMENTS.get(elementID)
            UTIL.DEBUG_OFF(ELEMENT)
            """add sectionID and elementID"""
            ELEMENT['ID']  = elementID 
            """INSTANCIATE ELM._Node objects"""
            instance = instanciate_element({elementID:ELEMENT}) 
            if instance == None: continue
            UTIL.DEBUG_OFF(instance)
            if isinstance(instance, ELM.Node):
                instances.append(instance)
        UTIL.DEBUG_OFF(instances)
        for instance in instances:
            lattice.add_node(instance)
        return lattice   # the complete lattice
    """ factory body -------- factory body -------- factory body -------- factory body -------- factory body -------- factory body -------- """
    with open(input_file,'r') as fileobject:
        try:
            in_data = yaml.load(fileobject,Loader=yaml.Loader)
        except Exception as ex:
            warnings.showwarning(
                    'InputError: {} - STOP'.format(str(ex)),
                    UserWarning,
                    'lattice_generator.py',
                    'factory()',
                    )
            sys.exit(1)
    fileobject.close()
    UTIL.DEBUG_OFF(in_data)

    # call lattice parser, get results
    results = doInputParser(in_data)

    flags = proces_flags(results.FLAGS)
    for k,v in flags.items(): UTIL.FLAGS[k] = v
    UTIL.DEBUG_OFF('global FLAGS after proces_flags():')
    UTIL.DEBUG_OFF(UTIL.FLAGS)

    parameters = proces_parameters(results.PARAMETERS)
    parameters['input_file'] = input_file
    for k,v in parameters.items(): UTIL.PARAMS[k] = v
    UTIL.DEBUG_OFF('global PARAMS after proces_parameters():')
    UTIL.DEBUG_OFF(UTIL.PARAMS)

    elements = proces_elements(results.ELEMENTS)
    UTIL.ELEMENTS = elements
    UTIL.DEBUG_OFF('ELEMENTS after proces_elements():')
    UTIL.DEBUG_OFF(UTIL.ELEMENTS)

    lat_elementIDs = results.LAT_ELMIDs
    UTIL.DEBUG_OFF('LAT_ELMIDs after proces_elements():')
    UTIL.DEBUG_OFF(lat_elementIDs)
    lattice = make_lattice(lat_elementIDs)
    # end of factory. return full Lattice object.
    return lattice
    
class TestLatticeGeneratorMethods(unittest.TestCase):
    def test_Lattice_Parser(self):
        print("\b----------------------------------------test_Lattice_Parser")
        input_file = "unittests/I5O200_19082022.yml"
        fileobject = open(input_file,'r')
        wfl = yaml.load(fileobject, Loader=yaml.FullLoader)
        fileobject.close()
        print('======================= yaml.dump(wfl,default_flow_style=True)')
        print(yaml.dump(wfl,default_flow_style=True))
        print('\n======================= wfl.items()')
        for i,v in iter(wfl.items()):
            print(i,' ==> ',v)
        
        seg = wfl['SEGMENTS']
        print("\n======================= seg = wfl['SEGMENTS']")
        print(seg)
        print('======================= segments')
        for i in seg:
            print(i)
        
        cell = wfl['CELLS']
        print("\n======================= seg = wfl['CELLS']")
        print(cell)
        print('======================= cells')
        for i in cell:
            print(i)
        
        sec = wfl['SECTIONS']
        print("\n======================= seg = wfl['SECTIONS']")
        print(sec)
        print('======================= sections')
        for i in sec:
            print(i)
        
        lattice = wfl['LATTICE']
        print('\n======================= lattice')
        for l in lattice:
            print(l)
    def test_Lattice_factory(self):
        print("\b----------------------------------------test_Lattice_factory")
        input_file = "unittests/I5O200_19082022.yml"
        lattice = factory(input_file)
        print(F'\u26dd  FINAL kinetic energy {lattice.seq[-1].ref_track[UTIL.EKOO]:.3f} [MeV] \u26dd')
        UTIL.DEBUG_ON('lattice_generator'   ,lattice.toString())
        UTIL.DEBUG_ON('SUMMARY in factory()',UTIL.SUMMARY)
        UTIL.DEBUG_ON('FLAGS in factory()'  ,UTIL.FLAGS)
        UTIL.DEBUG_ON('PARAMS in factory()' ,UTIL.PARAMS)

if __name__ == '__main__':
    unittest.main()
