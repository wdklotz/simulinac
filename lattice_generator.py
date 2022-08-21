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

def PRINT_PRETTY(obj):
    file = inspect.stack()[0].filename
    print(F'DEBUG_ON[{file}] ==> ',end="")
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON  = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

import setutil as util
import elements as ELM
import OXAL as OXA
import TTFG as TTF
import DYNG as DYN
from lattice import Lattice
from Ez0 import SFdata
from lattice_parser2 import parse as doInputParser
import PsMarkerAgent as psmkr

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
                instance = OXA.OXAL_G(ID,EzAvg,phiSoll,gap,freq,SFdata=util.PARAMS[fname],aperture=aperture,dWf=dWf)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            elif mapping == 'ttf':
                EzAvg = ELEMENT['EzAvg'] = util.PARAMS[fname].EzAvg
                instance = TTF.TTF_G(ID,EzAvg,phiSoll,gap,freq,SFdata=util.PARAMS[fname],aperture=aperture,dWf=dWf)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            elif mapping == 'dyn':
                EzAvg = ELEMENT['EzAvg'] = util.PARAMS[fname].EzAvg
                instance = DYN.DYN_G(ID,EzAvg,phiSoll,gap,freq,SFdata=util.PARAMS[fname],aperture=aperture,dWf=dWf)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
            else:
                instance = ELM.RFG(ID,EzAvg,phiSoll,gap,freq,mapping=mapping,aperture=aperture,dWf=dWf)
                ELEMENT['sec']    = attributes.get('sec','?')
                ELEMENT['EzPeak'] = EzPeak
        elif type == 'RFC':
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
            action = get_mandatory(attributes,'action',ID)
            if 'pspace' == action:
                which        = attributes.get('which','transvers')
                agent        = psmkr.PsMarkerAgent(which_action=which)
                instance     = ELM.MRK(ID,agents=[agent])
                instance.sec = attributes.get('sec','?')
                agent.set_parent(instance)
                DEBUG_OFF(instance.__dict__)
                DEBUG_OFF(agent.__dict__)

            # TODO !!
            # elif 'poincare' == action:
            #     prefix    = attributes['prefix']   if 'prefix'   in attributes else ''
            #     abszisse  = attributes['abscissa'] if 'abscissa' in attributes else 'z'
            #     ordinate  = attributes['ordinate'] if 'ordinate' in attributes else 'zp'
            #     instance = MRK.PoincareAction(label=label, prefix=prefix, abszisse=abszisse, ordinate=ordinate)
            #     instance['prefix']     = prefix
            #     instance['abszisse']   = abszisse
            #     instance['ordinate']   = ordinate
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
        if 'acccON'      in flags: util.FLAGS['accON']    = flags['accON']
        if 'periodic'    in flags: util.FLAGS['periodic'] = flags['periodic']
        if 'egf'         in flags: util.FLAGS['egf']      = flags['egf']
        if 'sigma'       in flags: util.FLAGS['sigma']    = flags['sigma']
        if 'KVout'       in flags: util.FLAGS['KVout']    = flags['KVout']
        if 'verbose'     in flags: util.FLAGS['verbose']  = flags['verbose']
        if 'useaper'     in flags: util.FLAGS['useaper']  = flags['useaper']
        if 'bucket'      in flags: util.FLAGS['bucket']   = flags['bucket']
        if 'csTrak'      in flags: util.FLAGS['csTrak']   = flags['csTrak']
        if 'marker'      in flags: util.FLAGS['marker']   = flags['marker']
        if 'pspace'      in flags: util.FLAGS['pspace']   = flags['pspace']
        if 'envelope'    in flags: util.FLAGS['envelope'] = flags['envelope']
        util.SUMMARY['accON'] = util.FLAGS.get('accON')
        if not util.FLAGS.get('accON'): util.FLAGS['dWf'] = 0.
        util.FLAGS['non_linear_mapping'] = False
        return flags
    def proces_parameters(parameters):
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
