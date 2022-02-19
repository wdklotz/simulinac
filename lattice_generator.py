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

def PRINT_PRETTY(obj):
    file = inspect.stack()[0].filename
    print('DEBUG_ON ==============>  '+file)
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON  = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

import setutil as util
import elements as ELM
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
    """ item: {ID:{attrinutes}} for each node """
    # def EzPeakToAverage(Ezpeak):
    #     # return 0.40 * EzPeak    # Pi mal Daumen: about 0.4*EzPeak from Superfish
    #     return EzPeak

    DEBUG_OFF(item)
    for ID,attributes in item.items():
        DEBUG_OFF(ID)
        DEBUG_OFF(attributes)
        # aperture   = PARAMS['aperture']    # default aperture
        type = attributes.get('type')
        if type == 'D':
            length       = get_mandatory(attributes,'length',ID)
            aperture     = attributes.get('aperture')
            instance     =  ELM.D(ID,length=length,aperture=aperture)
            instance.sec = attributes.get('sec','?')

        elif type == 'SIXD':
            length       = get_mandatory(attributes,'length',ID)
            aperture     = attributes.get('aperture')
            instance     = ELM.SIXD(ID,length=length,aperture=aperture)
            instance.sec = attributes.get('sec','?')

        elif type == 'QF':
            length        = get_mandatory(attributes,'length',ID)
            dBdz          = get_mandatory(attributes,"B'",ID)
            aperture      = get_mandatory(attributes,'aperture',ID)
            instance      = ELM.QF(ID,dBdz,length=length,aperture=aperture)
            instance.Bpole = dBdz*aperture      # Bpole
            instance.sec   = attributes.get('sec','?')

        elif type == 'QD':
            length         = get_mandatory(attributes,'length',ID)
            dBdz           = get_mandatory(attributes,"B'",ID)
            aperture       = get_mandatory(attributes,'aperture',ID)
            instance       = ELM.QD(ID,dBdz,length=length,aperture=aperture)
            instance.Bpole = dBdz*aperture      # Bpole
            instance.sec   = attributes.get('sec','?')

        elif type == 'RFG':
            phiSoll   = radians(get_mandatory(attributes,"PhiSync",ID))
            freq      = float(get_mandatory(attributes,"freq",ID))
            gap       = get_mandatory(attributes,'gap',ID)
            aperture  = get_mandatory(attributes,'aperture',ID)
            dWf       = util.FLAGS['dWf']
            mapping   = get_mandatory(attributes,'mapping',ID)
            EzPeak    = get_mandatory(attributes,"EzPeak",ID)
            mapping   = attributes.get('mapping','t3d')
            if mapping == 'ttf' or mapping == 'dyn' or mapping == 'oxal': # SF-data
                fname = get_mandatory(attributes,"SFdata",ID)
                if fname not in util.PARAMS:
                    gap_cm = gap*100     # Watch out!
                    util.PARAMS[fname] = SFdata(fname,EzPeak=EzPeak,gap=gap_cm)
                EzAvg    = util.PARAMS[fname].EzAvg
                instance = ELM.RFG(ID,EzAvg,phiSoll,gap,freq,mapping=mapping,SFdata=util.PARAMS[fname],aperture=aperture,dWf=dWf)
            else:
                EzAvg    = EzPeak
                instance = ELM.RFG(ID,EzAvg,phiSoll,gap,freq,mapping=mapping,aperture=aperture,dWf=dWf)
            element           = util.ELEMENTS[ID]
            element['SFdata'] = instance.SFdata
            element['EzAvg']  = EzAvg
            instance.sec      = attributes.get('sec','?')

        elif type == 'RFC':
            label     = attributes['ID']
            PhiSoll   = radians(get_mandatory(attributes,"PhiSync",label))
            freq      = float(get_mandatory(attributes,"freq",label))
            gap       = get_mandatory(attributes,'gap',label)
            aperture  = get_mandatory(attributes,'aperture',label)
            dWf       = util.FLAGS['dWf']
            length    = get_mandatory(attributes,'length',label)
            mapping   = get_mandatory(attributes,'mapping',label)
            EzPeak    = get_mandatory(attributes,"EzPeak",label)
            # EzAvg     = attributes['EzAvg'] if 'EzAvg' in attributes else EzPeakToAverage(EzPeak)
            EzAvg     = attributes.get('EzAvg',EzPeakToAverage(EzPeak))
            if mapping == None:
                mapping = 't3d'
            if mapping == 'ttf' or mapping == 'dyn' or mapping == 'oxal': # SF-data
                fname     = get_mandatory(attributes,"SFdata",label)
                if fname not in util.PARAMS:
                    gap_cm = gap*100     # Watch out!
                    util.PARAMS[fname] = SFdata(fname,EzPeak=EzPeak,gap=gap_cm)
                EzAvg = util.PARAMS[fname].EzAvg
                instance  =  ELM.RFC(EzAvg=EzAvg,label=label,PhiSoll=PhiSoll,fRF=freq,gap=gap,aperture=aperture,dWf=dWf,length=length,mapping=mapping,SFdata=util.PARAMS[fname])
                pass
            else:
                # EzAvg = EzPeakToAverage(EzPeak)  TODO ??
                instance  =  ELM.RFC(EzAvg=EzAvg,label=label,PhiSoll=PhiSoll,fRF=freq,gap=gap,aperture=aperture,dWf=dWf,length=length,mapping=mapping)
            element = util.ELEMENTS[label]
            element['EzAvg']     = EzAvg
            instance.sec = attributes.get('sec','?')
            instance['EzAvg']    = EzAvg
            instance['EzPeak']   = EzPeak
            instance['label']    = label
            instance['PhiSoll']  = PhiSoll
            instance['freq']     = freq
            instance['gap']      = gap
            instance['aperture'] = aperture
            instance['dWf']      = dWf
            instance['length']   = length
            instance['mapping']  = mapping

        elif type == 'GAP':
            gap       = get_mandatory(attributes,'gap',ID)
            EzAvg     = get_mandatory(attributes,"EzAvg",ID)
            phiSoll   = radians(get_mandatory(attributes,"PhiSync",ID))
            freq      = float(get_mandatory(attributes,"freq",ID))
            dWf       = util.FLAGS['dWf']
            aperture  = get_mandatory(attributes,'aperture',ID)
            instance  =  ELM.GAP(ID,EzAvg,phiSoll,gap,freq,aperture=aperture,dWf=dWf)
            element = util.ELEMENTS[ID]
            element['EzPeak'] = EzAvg
            instance.EzPeak   = EzPeak
            instance.sec      = attributes.get('sec','?')

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
        if 'accON' in flags:
            if flags['accON']:
                util.FLAGS['dWf'] = 1.
                util.SUMMARY['accON'] = True
            else:
                util.FLAGS['dWf'] = 0.
                util.SUMMARY['accON'] = False
        if 'periodic'    in flags: util.FLAGS['periodic'] = flags['periodic']
        if 'egf'         in flags: util.FLAGS['egf']      = flags['egf']
        if 'sigma'       in flags: util.FLAGS['sigma']    = flags['sigma']
        if 'KVout'       in flags: util.FLAGS['KVout']    = flags['KVout']
        if 'verbose'     in flags: util.FLAGS['verbose']  = flags['verbose']
        if 'express'     in flags: util.FLAGS['express']  = flags['express']
        if 'useaper'     in flags: util.FLAGS['useaper']  = flags['useaper']
        if 'bucket'      in flags: util.FLAGS['bucket']   = flags['bucket']
        if 'csTrak'      in flags: util.FLAGS['csTrak']   = flags['csTrak']
        if 'marker'      in flags: util.FLAGS['marker']   = flags['marker']
        if 'pspace'      in flags: util.FLAGS['pspace']   = flags['pspace']
        return flags
    def proces_parameters(parameters):
        """ fills global PARAMETERS"""
        if 'Tkin'             in parameters: util.PARAMS['injection_energy'] = parameters['Tkin']
        if 'phi_sync'         in parameters: util.PARAMS['phisoll']          = parameters['phi_sync']
        if 'gap'              in parameters: util.PARAMS['gap']              = parameters['gap']
        if 'cav_len'          in parameters: util.PARAMS['cavity_laenge']    = parameters['cav_len']
        if 'frequency'        in parameters: util.PARAMS['frequency']        = parameters['frequency']
        if 'ql'               in parameters: util.PARAMS['ql']               = parameters['ql']
        if 'windings'         in parameters: util.PARAMS['nbwindgs']         = parameters['windings']
        if 'nbsigma'          in parameters: util.PARAMS['nbsigma']          = parameters['nbsigma']
        if 'aperture'         in parameters: util.PARAMS['aperture']         = parameters['aperture'] 
        if 'emitx_i'          in parameters: util.PARAMS['emitx_i']          = parameters['emitx_i']
        if 'emity_i'          in parameters: util.PARAMS['emity_i']          = parameters['emity_i']
        if 'betax_i'          in parameters: util.PARAMS['betax_i']          = parameters['betax_i']
        if 'betay_i'          in parameters: util.PARAMS['betay_i']          = parameters['betay_i']
        if 'alfax_i'          in parameters: util.PARAMS['alfax_i']          = parameters['alfax_i']
        if 'alfay_i'          in parameters: util.PARAMS['alfay_i']          = parameters['alfay_i']
        if 'mapping'          in parameters: util.PARAMS['mapping']          = parameters['mapping']
        if 'thins'            in parameters: util.PARAMS['thins']            = parameters['thins']
        if 'DT2T'             in parameters: util.PARAMS['DT2T']             = parameters['DT2T']
        if 'lattvers'         in parameters: util.PARAMS['lattice_version']  = parameters['lattvers']
        return parameters
    def proces_elements(elements):
        """fills global ELEMENTS"""
        util.ELEMENTS = elements
        return elements
    def make_lattice(elementIDs):
        # DEBUG_ON(elementIDs)
        lattice = Lattice()
        instances = []
        for elementID in elementIDs:
            # print("A")
            # DEBUG_ON(elementID)
            element = util.ELEMENTS.get(elementID)
            # print("B")
            # DEBUG_ON(element)
            """add sectionID and elementID"""
            element['ID']  = elementID 
            # repack {ID:{attributes}} for instanciate_element(...)
            item = {elementID:element}
            """INSTANCIATE ELM._Node objects"""
            instance = instanciate_element(item)
            # print("C")
            # DEBUG_ON(instance)
            if isinstance(instance, ELM.Node):
                # lattice.add_element(instance)
                instances.append(instance)
            # list of thin quad instances
            elif isinstance(instance,list):
                # [lattice.add_element(x) for x in instance]
                instances += instance
        # print("D")
        # DEBUG_ON(instances)
        for instance in instances:
            # print('E')
            # DEBUG_ON(instance.__dict__)
            lattice.add_element(instance)
        return lattice   # the complete lattice
    ## factory body -------- factory body -------- factory body -------- factory body -------- factory body -------- factory body --------
    ## factory body -------- factory body -------- factory body -------- factory body -------- factory body -------- factory body --------
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

    results = doInputParser(in_data)  # call lattice parser, get results

    flags = proces_flags(results.FLAGS)
    DEBUG_OFF('global FLAGS after proces_flags():')
    DEBUG_OFF(util.FLAGS)

    parameters = proces_parameters(results.PARAMETERS)
    DEBUG_OFF('global PARAMS after proces_parameters():')
    DEBUG_OFF(util.PARAMS)
    util.PARAMS['sollteilchen'](tkin=util.PARAMS['injection_energy'])
    
    elements = proces_elements(results.ELEMENTS)
    DEBUG_OFF('ELEMENTS after proces_elements():')
    DEBUG_OFF(util.ELEMENTS)

    lat_elementIDs = results.LAT_ELMIDs
    lattice = make_lattice(lat_elementIDs)
    DEBUG_OFF('lattice_generator >>{}'.format(lattice.toString()))
    util.SUMMARY['lattice length [m]'] = util.PARAMS['lattice_length']  = lattice.length
    DEBUG_OFF('SUMMARY in factory() {}'.format(util.SUMMARY))
    # end of factory(...)
    return lattice
def test0(input_file):
    print('---------------------------------TEST0')
    wfl= []
    fileobject = open(input_file,'r')
    wfl = yaml.load(fileobject)
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
# def test1(input_file):
#     print('---------------------------------TEST1')
#     lattice = factory(input_file)
#     print('%%%%%%%%%%%%%%%%%%% Left------->Right')
#     for cnt,node in enumerate(iter(lattice)):
#         print(cnt,'{:38s} {:38s}'.format(repr(node),repr(node.next)))
#     print('\n%%%%%%%%%%%%%%%%%%% Right------->Left')
#     lattice.toggle_iteration()
#     for cnt,node in enumerate(iter(lattice)):
#         print(cnt,'{:38s} {:38s}'.format(repr(node),repr(node.prev)))

if __name__ == '__main__':
    input_file = "REF-Läufe\DG6FG6D-v10.0.1-ref\simuIN.yml"
    # input_file = 'yml/v10.0.0_compat_IN.yml'
    # input_file = 'yml/simuIN.yml'
    # input_file = 'nwlat/nwlatIN.yml'
    test0(input_file)
    # test1(input_file)
