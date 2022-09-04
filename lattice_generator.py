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
import sys
from math import radians
import yaml
import warnings
import pprint, inspect
import unittest

import setutil as util
import elements as ELM
import OXAL as OXA
import TTFG as TTF
import DYNG as DYN
import PsMarkerAgent as psmkr
import PoincareMarkerAgent as pcmkr
from lattice import Lattice
from Ez0 import SFdata
from lattice_parser2 import parse as doInputParser

def PRINT_PRETTY(obj):
    file = inspect.stack()[0].filename
    print(F'DEBUG_ON[{file}] ==> ',end="")
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON  = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

def make_counter():
    count = 0
    def inner():
        nonlocal count
        count += 1
        return count
    return inner
counter1 = make_counter()
counter2 = make_counter()

def check_marker_incompatible_with(prog,ID):
    ret = False
    this_prog = sys.argv[0]
    if prog != this_prog:
        print(util.colors.RED+f'WARN: Marker {ID} incompatible with {this_prog}. Will be skipped'+util.colors.ENDC)
        ret = True
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
        DEBUG_OFF(F"ID={ID} attributes={attributes}")
        ELEMENT = util.ELEMENTS[ID]          # the item in the ELEMENT list
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
            dWf       = util.FLAGS['dWf']
            EzPeak    = get_mandatory(attributes,"EzPeak",ID)
            mapping   = attributes.get('mapping','t3d')
            if mapping == 'ttf' or mapping == 'dyn' or mapping == 'oxal': # SF-data
                fname = get_mandatory(attributes,"SFdata",ID)
                if fname not in util.PARAMS:
                    gap_cm = gap*100     # Watch out!
                    util.PARAMS[fname] = SFdata(fname,EzPeak=EzPeak,gap=gap_cm)
            else:
                ELEMENT['EzAvg']  = EzAvg = EzPeak
                ELEMENT['SFdata'] = None
            if mapping == 'oxal':
                EzAvg = ELEMENT['EzAvg'] = util.PARAMS[fname].EzAvg
                instance = OXA.OXAL_G(ID,EzAvg,phiSoll,gap,freq,SFdata=util.PARAMS[fname],particle=util.Proton(util.PARAMS['injection_energy']),position=(0.,0.,0.),aperture=aperture,dWf=dWf)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            elif mapping == 'ttf':
                EzAvg = ELEMENT['EzAvg'] = util.PARAMS[fname].EzAvg
                instance = TTF.TTF_G(ID,EzAvg,phiSoll,gap,freq,SFdata=util.PARAMS[fname],particle=util.Proton(util.PARAMS['injection_energy']),position=(0.,0.,0.),aperture=aperture,dWf=dWf)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            elif mapping == 'dyn':
                if counter1() < 5:
                    print(util.colors.RED+"WARN: dyn mapping is broken"+util.colors.ENDC)
                else:
                    pass
                EzAvg = ELEMENT['EzAvg'] = util.PARAMS[fname].EzAvg
                instance = DYN.DYN_G(ID,EzAvg,phiSoll,gap,freq,SFdata=util.PARAMS[fname],particle=util.Proton(util.PARAMS['injection_energy']),position=(0.,0.,0.),aperture=aperture,dWf=dWf)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            else:
                instance = ELM.RFG(ID,EzAvg,phiSoll,gap,freq,SFdata=None,particle=util.Proton(util.PARAMS['injection_energy']),position=(0.,0.,0.),aperture=aperture,dWf=dWf,mapping=mapping)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
        elif type == 'RFC':
            if counter2() < 5:
                print(util.colors.RED+"WARN: RFC node is broken"+util.colors.ENDC)
            else:
                sys.exit(1)
            phiSoll   = radians(get_mandatory(attributes,"PhiSync",ID))
            freq      = float(get_mandatory(attributes,"freq",ID))
            gap       = get_mandatory(attributes,'gap',ID)
            length    = get_mandatory(attributes,'length',ID)
            aperture  = get_mandatory(attributes,'aperture',ID)
            dWf       = util.FLAGS['dWf']
            EzPeak    = get_mandatory(attributes,"EzPeak",ID)
            mapping   = attributes.get('mapping','t3d')
            if mapping == 'ttf' or mapping == 'dyn' or mapping == 'oxal': # SF-data
                fname = get_mandatory(attributes,"SFdata",ID)
                if fname not in util.PARAMS:
                    gap_cm = gap*100     # Watch out!
                    util.PARAMS[fname] = SFdata(fname,EzPeak=EzPeak,gap=gap_cm)
            else:
                ELEMENT['EzAvg']  = EzAvg = EzPeak
                ELEMENT['SFdata'] = None
            if mapping == 'oxal':
                EzAvg = ELEMENT['EzAvg'] = util.PARAMS[fname].EzAvg
                instance = OXA.OXAL_C(ID,EzAvg,phiSoll,gap,length,freq,SFdata=util.PARAMS[fname],aperture=aperture,dWf=dWf)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            elif mapping == 'ttf':
                EzAvg = ELEMENT['EzAvg'] = util.PARAMS[fname].EzAvg
                instance = TTF.TTF_C(ID,EzAvg,phiSoll,gap,length,freq,SFdata=util.PARAMS[fname],aperture=aperture,dWf=dWf)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            elif mapping == 'dyn':
                EzAvg = ELEMENT['EzAvg'] = util.PARAMS[fname].EzAvg
                instance = DYN.DYN_C(ID,EzAvg,phiSoll,gap,length,freq,SFdata=util.PARAMS[fname],aperture=aperture,dWf=dWf)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            else:
                instance = ELM.RFC(ID,EzAvg,phiSoll,gap,freq,length,mapping=mapping,aperture=aperture,dWf=dWf)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
        elif type == 'GAP':
            gap       = get_mandatory(attributes,'gap',ID)
            EzPeak    = get_mandatory(attributes,"EzPeak",ID)
            phiSoll   = radians(get_mandatory(attributes,"PhiSync",ID))
            freq      = float(get_mandatory(attributes,"freq",ID))
            dWf       = util.FLAGS['dWf']
            aperture  = get_mandatory(attributes,'aperture',ID)
            EzAvg     = EzPeak
            instance  =  ELM.GAP(ID,EzAvg,phiSoll,gap,freq,aperture=aperture,dWf=dWf)
            ELEMENT['EzPeak'] = EzPeak
            ELEMENT['sec']    = attributes.get('sec','?')
            ELEMENT['EzAvg']  = EzAvg
        elif type == 'MRK':
            active  = attributes.get('active',False)
            if not active: continue       # exclude this marker from lattice
            ELEMENT = util.ELEMENTS[ID]
            action  = get_mandatory(attributes,'action',ID)
            if action == 'pspace':
                # A marker for simu.py ?
                if check_marker_incompatible_with('simu.py',ID): continue    # exclude this marker from lattice
                agent    = psmkr.PsMarkerAgent()
                instance = ELM.MRK(ID,agent,active)
                agent.set_parent(instance)
                sec = attributes.get('sec','?') 
                ELEMENT['sec']   = sec
                DEBUG_OFF(ELEMENT)
                DEBUG_OFF(instance.toString())
                DEBUG_OFF(agent.__dict__)

            elif action == 'pcrcut':
                # A marker for tracker.py ?
                if check_marker_incompatible_with('tracker.py',ID): continue    # exclude this marker from lattice
                sec        = attributes.get('sec','?') 
                prefix     = attributes.get('prefix','frames')
                abscissa   = attributes.get('abscissa','z')
                ordinate   = attributes.get('ordinate','zp')
                agent      = pcmkr.PoincareMarkerAgent(sec,prefix,abscissa,ordinate)
                instance   = ELM.MRK(ID,agent,active)
                agent.set_parent(instance)
                ELEMENT['sec']      = sec
                ELEMENT['prefix']   = prefix
                ELEMENT['abscissa'] = abscissa
                ELEMENT['ordinate'] = ordinate
                DEBUG_OFF(ELEMENT)
                DEBUG_OFF(instance.__dict__)
                DEBUG_OFF(agent.__dict__)
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
    """ factory creates a lattice from input-file """

    def proces_flags(flags):
        """fills global FLAGS"""        
        util.FLAGS['accON']    = flags.get('accON',True)               # acceleration ON
        util.FLAGS['periodic'] = flags.get('periodic',False)           # periodic lattice? default
        util.FLAGS['egf']      = flags.get('egf',False)                # emittance grow flag default
        util.FLAGS['sigma']    = flags.get('sigma',True)               # beam sizes by sigma-tracking
        util.FLAGS['KVout']    = flags.get('KVout',False)              # print a dictionary of Key-Value pairs, no display
        util.FLAGS['verbose']  = flags.get('verbose',0)                # print flag default = 0
        util.FLAGS['useaper']  = flags.get('useaper',False)            # use aperture check for quads and rf-gaps
        util.FLAGS['bucket']   = flags.get('bucket',False)             # plot bucket
        util.FLAGS['csTrak']   = flags.get('csTrak',True)              # plot CS trajectories
        util.FLAGS['maction']  = flags.get('maction',False)            # call marker actions
        util.FLAGS['pspace']   = flags.get('pspace',False)             # plot CS twiss ellipses at entrance
        util.FLAGS['envelope'] = flags.get('envelope',False)           # plot transverse envelopes
        util.FLAGS['dWf']      = 1 if util.FLAGS.get('accON') else 0   # acceleration on/off flag 1=on,0=off
        util.SUMMARY['accON']  = util.FLAGS.get('accON')
        util.FLAGS['non_linear_mapping'] = False
        return flags
    def proces_parameters(parameters):   #TODO use dict.get()
        """ fills global PARAMETERS"""
        if 'Tkin'             in parameters: util.PARAMS['injection_energy'] = parameters['Tkin']
        if 'DT2T'             in parameters: util.PARAMS['DT2T']             = parameters['DT2T']
        if 'emitw'            in parameters: util.PARAMS['emitw']            = parameters['emitw']
        if 'Dphi0'            in parameters: util.PARAMS['Dphi0']            = parameters['Dphi0']
        if 'emitx_i'          in parameters: util.PARAMS['emitx_i']          = parameters['emitx_i']
        if 'emity_i'          in parameters: util.PARAMS['emity_i']          = parameters['emity_i']
        if 'betax_i'          in parameters: util.PARAMS['betax_i']          = parameters['betax_i']
        if 'betay_i'          in parameters: util.PARAMS['betay_i']          = parameters['betay_i']
        if 'alfax_i'          in parameters: util.PARAMS['alfax_i']          = parameters['alfax_i']
        if 'alfay_i'          in parameters: util.PARAMS['alfay_i']          = parameters['alfay_i']
        # if 'alfaw_i'          in parameters: util.PARAMS['alfaw_i']          = parameters['alfaw_i']
        if 'nbsigma'          in parameters: util.PARAMS['nbsigma']          = parameters['nbsigma']
        if 'lattvers'         in parameters: util.PARAMS['lattice_version']  = parameters['lattvers']
        if 'mapping'          in parameters: util.PARAMS['mapping']          = parameters['mapping']
        """
        if 'frequency'        in parameters: util.PARAMS['frequency']        = parameters['frequency']
        if 'phi_sync'         in parameters: util.PARAMS['phisoll']          = parameters['phi_sync']
        if 'gap'              in parameters: util.PARAMS['gap']              = parameters['gap']
        if 'cav_len'          in parameters: util.PARAMS['cavity_laenge']    = parameters['cav_len']
        if 'ql'               in parameters: util.PARAMS['ql']               = parameters['ql']
        if 'windings'         in parameters: util.PARAMS['nbwindgs']         = parameters['windings']
        if 'aperture'         in parameters: util.PARAMS['aperture']         = parameters['aperture'] 
        if 'thins'            in parameters: util.PARAMS['thins']            = parameters['thins']
        """
        return parameters
    def proces_elements(elements):
        """fills global ELEMENTS"""
        util.ELEMENTS = elements
        return elements
    def make_lattice(elementIDs,injection_energy):
        # DEBUG_ON(elementIDs)
        lattice = Lattice(injection_energy)
        instances = []
        for elementID in elementIDs:
            # print("A"); DEBUG_ON(elementID)
            ELEMENT = util.ELEMENTS.get(elementID)
            if ELEMENT.get('mapping',"") in ['base','ttf','dyn']:   # non_linear_mapping in lattice?
                util.FLAGS['non_linear_mapping'] = True
            # print("B"); DEBUG_ON(element)
            """add sectionID and elementID"""
            ELEMENT['ID']  = elementID 
            # repack {ID:{attributes}} for instanciate_element(...)
            item = {elementID:ELEMENT}
            """INSTANCIATE ELM._Node objects"""
            instance = instanciate_element(item) 
            if instance == None: continue
            # print("C"); DEBUG_ON(instance)
            if isinstance(instance, ELM.Node):
                # lattice.add_node(instance)
                instances.append(instance)
            # list of thin quad instances
            elif isinstance(instance,list):
                # [lattice.add_node(x) for x in instance]
                instances += instance
        # print("D"); DEBUG_ON(instances)
        for instance in instances:
            lattice.add_node(instance)
        return lattice   # the complete lattice


    ## factory body -------- factory body -------- factory body -------- factory body -------- factory body -------- factory body --------
    util.SUMMARY['input file'] = util.PARAMS['input_file'] = input_file
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
    DEBUG_OFF(in_data)

    # call lattice parser, get results
    results = doInputParser(in_data)

    flags = proces_flags(results.FLAGS)
    DEBUG_OFF('global FLAGS after proces_flags():')
    DEBUG_OFF(util.FLAGS)

    parameters = proces_parameters(results.PARAMETERS)
    DEBUG_OFF('global PARAMS after proces_parameters():')
    DEBUG_OFF(util.PARAMS)
    
    elements = proces_elements(results.ELEMENTS)
    DEBUG_OFF('ELEMENTS after proces_elements():')
    DEBUG_OFF(util.ELEMENTS)

    lat_elementIDs = results.LAT_ELMIDs
    lattice = make_lattice(lat_elementIDs,util.PARAMS['injection_energy'])
    DEBUG_OFF('lattice_generator >>{}'.format(lattice.toString()))
    DEBUG_OFF('SUMMARY in factory() {}'.format(util.SUMMARY))
    # end of factory(...)
    return lattice
    
class TestLatticeGeneratorMethods(unittest.TestCase):
    def test_Lattice_Parser(self):
        print("\b----------------------------------------test_Lattice_Parser")
        input_file = "REF-Läufe\DG6FG6D-v10.0.1-ref\simuIN.yml"
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
        print('\n======================= segment')
        for i in seg:
            print(i)
        lattice = wfl['LATTICE']
        print('\n======================= lattice')
        for l in lattice:
            print(l)

if __name__ == '__main__':
    unittest.main()
