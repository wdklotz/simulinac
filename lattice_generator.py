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

import setutil as util
import elements as ELM
from lattice import Lattice
from Ez0 import SFdata
import marker_actions as MRK
from new_lattice_parser import parse
from collections import namedtuple

DEBUG_ON  = util.DEB.get('ON')
DEBUG_OFF = util.DEB.get('OFF')

# parse and generate latttice
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

#---- recursive version
# def flatten(lis):
#     """Given a list, possibly nested to any level, return it flattened."""
#     new_lis = []
#     for item in lis:
#         if type(item) == type([]):
#             new_lis.extend(flatten(item))
#         else:
#             new_lis.append(item)
#     return new_lis

# def liofd2d(l):
#     """Transform list of dicts to dict"""
#     return {k:v for d in l for k,v in d.items()}

def replace_QF_with_QFth_lattice(slices,k0,length,label,particle,aperture):
    lattice = Lattice()
    thinlen = length/slices
    for nb in range(slices):
        if util.FLAGS['express']:
            instance = ELM.QFthx(k0=k0,length=thinlen,label=label,particle=particle,aperture=aperture)
        else:
            instance = ELM.QFth(k0=k0,length=thinlen,label=label,particle=particle,aperture=aperture)
        lattice.add_element(instance)
    return lattice

def replace_QD_with_QDth_lattice(slices,k0,length,label,particle,aperture):
    lattice = Lattice()
    thinlen = length/slices
    for nb in range(slices):
        if util.FLAGS['express']:
            instance = ELM.QDthx(k0=k0,length=thinlen,label=label,particle=particle,aperture=aperture)
        else:
            instance = ELM.QDth(k0=k0,length=thinlen,label=label,particle=particle,aperture=aperture)
        lattice.add_element(instance)
    return lattice
def instanciate_element(item):
    """ item: {ID:{attrinutes}} for each element """
    def EzPeakToAverage(Ezpeak):
        return 0.78 * EzPeak    # ~0.748 * EzPeak from Superfish

    DEBUG_OFF(item)
    for ID,attributes in item.items():
        DEBUG_OFF(ID)
        DEBUG_OFF(attributes)
        # aperture   = PARAMS['aperture']    # default aperture
        key = attributes.get('type')
        if key == 'D':
            label    = ID
            length   = get_mandatory(attributes,'length',label)
            aperture = attributes['aperture'] if 'aperture' in attributes else None
            instance =  ELM.D(length=length,label=label,aperture=aperture)
            instance['label']    = label
            instance['length']   = length
            instance['aperture'] = aperture

        elif key == 'SIXD':
            label     = attributes['ID']+'#'
            length    = get_mandatory(attributes,'length',label)
            aperture = attributes['aperture'] if 'aperture' in attributes else None
            instance  = ELM.SIXD(length=length,label=label,aperture=aperture)
            instance['label']    = label
            instance['length']   = length
            instance['aperture'] = aperture

        elif key == 'QF':
            label    = attributes['ID']
            length   = get_mandatory(attributes,'length',label)
            dBdz     = get_mandatory(attributes,"B'",label)
            slices   = get_mandatory(attributes,'slices',label)
            aperture = get_mandatory(attributes,'aperture',label)
            kq       = dBdz/util.PARAMS['sollteilchen'].brho
            Bpole    = dBdz*aperture
            if slices > 1:
                instance = replace_QF_with_QFth_lattice(slices,kq,length,label,util.PARAMS['sollteilchen'],aperture)
            elif slices <= 1:
                instance = ELM.QF(k0=kq,length=length,label=label,aperture=aperture)
            instance['label']  = label
            instance['dBdz']   = dBdz
            instance['bore']   = aperture
            instance['Bpole']  = Bpole
            pass

        elif key == 'QD':
            label    = attributes['ID']
            length   = get_mandatory(attributes,'length',label)
            dBdz     = get_mandatory(attributes,"B'",label)
            slices   = get_mandatory(attributes,'slices',label)
            aperture = get_mandatory(attributes,'aperture',label)
            kq       = dBdz/util.PARAMS['sollteilchen'].brho
            Bpole    = dBdz*aperture
            if slices > 1:
                instance = replace_QD_with_QDth_lattice(slices,kq,length,label,util.PARAMS['sollteilchen'],aperture)
            elif slices <= 1:
                instance = ELM.QD(k0=kq,length=length,label=label,aperture=aperture)
            instance['label']  = label
            instance['dBdz']   = dBdz
            instance['bore']   = aperture
            instance['Bpole']  = Bpole

        elif key == 'RFG':
            label     = attributes['ID']
            PhiSoll   = radians(get_mandatory(attributes,"PhiSync",label))
            freq      = float(get_mandatory(attributes,"freq",label))
            gap       = get_mandatory(attributes,'gap',label)
            aperture  = get_mandatory(attributes,'aperture',label)
            dWf       = util.FLAGS['dWf']
            mapping   = get_mandatory(attributes,'mapping',label)
            EzPeak    = get_mandatory(attributes,"EzPeak",label)
            EzAvg     = attributes['EzAvg'] if 'EzAvg' in attributes else EzPeakToAverage(EzPeak)
            if mapping == None:
                mapping = 't3d'
                # EzAvg = EzPeakToAverage(EzPeak)
            if mapping == 'ttf' or mapping == 'dyn' or mapping == 'oxal': # SF-data
                fname     = get_mandatory(attributes,"SFdata",label)
                if fname not in util.PARAMS:
                    half_gap_cm = gap*50     # Watch out!
                    util.PARAMS[fname] = SFdata(fname,EzPeak=EzPeak,gap=half_gap_cm)
                EzAvg = util.PARAMS[fname].EzAvg
                instance  =  ELM.RFG(EzAvg=EzAvg,PhiSoll=PhiSoll,fRF=freq,label=label,gap=gap,mapping=mapping,dWf=dWf,aperture=aperture,SFdata=util.PARAMS[fname])
                pass
            else:
                # EzAvg = EzPeakToAverage(EzPeak)
                instance  = ELM.RFG(EzAvg=EzAvg,PhiSoll=PhiSoll,fRF=freq,label=label,gap=gap,mapping=mapping,dWf=dWf,aperture=aperture)
            element = util.ELEMENTS[label]
            element['EzAvg']     = EzAvg
            instance['EzAvg']    = EzAvg
            instance['EzPeak']   = EzPeak
            instance['label']    = label
            instance['PhiSoll']  = PhiSoll
            instance['freq']     = freq
            instance['gap']      = gap
            instance['aperture'] = aperture
            instance['dWf']      = dWf
            instance['mapping']  = mapping

        elif key == 'RFC':
            label     = attributes['ID']
            PhiSoll   = radians(get_mandatory(attributes,"PhiSync",label))
            freq      = float(get_mandatory(attributes,"freq",label))
            gap       = get_mandatory(attributes,'gap',label)
            aperture  = get_mandatory(attributes,'aperture',label)
            dWf       = util.FLAGS['dWf']
            length    = get_mandatory(attributes,'length',label)
            mapping   = get_mandatory(attributes,'mapping',label)
            EzPeak    = get_mandatory(attributes,"EzPeak",label)
            EzAvg     = attributes['EzAvg'] if 'EzAvg' in attributes else EzPeakToAverage(EzPeak)
            if mapping == None:
                mapping = 't3d'
                # EzAvg = EzPeakToAverage(EzPeak)
            if mapping == 'ttf' or mapping == 'dyn' or mapping == 'oxal': # SF-data
                fname     = get_mandatory(attributes,"SFdata",label)
                if fname not in util.PARAMS:
                    half_gap_cm = gap*50     # Watch out!
                    util.PARAMS[fname] = SFdata(fname,EzPeak=EzPeak,gap=half_gap_cm)
                EzAvg = util.PARAMS[fname].EzAvg
                instance  =  ELM.RFC(EzAvg=EzAvg,label=label,PhiSoll=PhiSoll,fRF=freq,gap=gap,aperture=aperture,dWf=dWf,length=length,mapping=mapping,SFdata=util.PARAMS[fname])
                pass
            else:
                # EzAvg = EzPeakToAverage(EzPeak)
                instance  =  ELM.RFC(EzAvg=EzAvg,label=label,PhiSoll=PhiSoll,fRF=freq,gap=gap,aperture=aperture,dWf=dWf,length=length,mapping=mapping)
            element = util.ELEMENTS[label]
            element['EzAvg']     = EzAvg
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

        elif key == 'GAP':
            label     = attributes['ID']
            gap       = get_mandatory(attributes,'gap',label)
            EzPeak    = get_mandatory(attributes,"EzPeak",label)
            EzAvg     = EzPeak
            PhiSoll   = radians(get_mandatory(attributes,"PhiSync",label))
            freq      = float(get_mandatory(attributes,"freq",label))
            dWf       = util.FLAGS['dWf']
            aperture  = get_mandatory(attributes,'aperture',label)
            instance  =  ELM.GAP(EzAvg=EzAvg,PhiSoll=PhiSoll,fRF=freq,label=label,gap=gap,dWf=dWf,aperture=aperture)
            instance['EzAvg']   = EzAvg
            instance['EzPeak']  = EzPeak
            instance['label']   = label
            instance['gap']     = gap
            instance['PhiSoll'] = PhiSoll
            instance['freq']    = freq
            instance['dWf']     = dWf

        elif key == 'MRK':
            label     = attributes['ID']
            action    = get_mandatory(attributes,'action',label)
            if 'poincare' == action:
                prefix    = attributes['prefix'] if 'prefix' in attributes else ''
                abszisse  = attributes['abscissa'] if 'abscissa' in attributes else 'z'
                ordinate  = attributes['ordinate'] if 'ordinate' in attributes else 'zp'
                instance = MRK.PoincareAction(label=label, prefix=prefix, abszisse=abszisse, ordinate=ordinate)
                instance['prefix']     = prefix
                instance['abszisse']   = abszisse
                instance['ordinate']   = ordinate
            else:
                instance  = ELM.MRK(label=label,action=action)
            instance['label']      = label
            instance['action']     = action
        else:
            warnings.showwarning(
                    'InputError: Unknown element type encountered: "{}" - STOP'.format(key),
                    UserWarning,
                    'lattice_generator.py',
                    'instanciate_element()',
                    )
            sys.exit(1)
        try:
            sec = attributes['sec']    #can fail because sections are not mandatory
        except:
            pass
        else:
            instance.section = sec
        break
    return instance


def factory_new(input_file):
    """ factory creates a lattice from input-file """
    #--------
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
        if 'pspace'      in flags: util.FLAGS['pspace']   = flags['pspace']
        return flags
    # --------
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
        if 'DT2T'             in parameters: util.PARAMS['DT2T']             = parameters['DT2T']
        if 'lattvers'         in parameters: util.PARAMS['lattice_version']  = parameters['lattvers']
        return parameters
    #--------
    def proces_elements(elements):
        """fills global ELEMENTS"""
        util.ELEMENTS = elements
        return elements
    # --------
    def proces_sections(results):
        """fills global SECTIONS"""
        sections = results.SECTIONSx
        util.SECTIONS={"allIDs":results.SECTIONSx, "uniqueIDs":results.SECTIONSu}
        return sections
    # --------
    def proces_lattice(lattice):
        """fills global LATTICE"""
        util.LATTICE = lattice
        return lattice
    # --------
    def make_lattice():
        lattice = Lattice()
        sections = util.SECTIONS['allIDs'] # dict of sections and their expanded key-list z.B. {LE:[id1,id2,...],He:[id1,id2,...],...}
        DEBUG_OFF(sections)
        sectionIDs = util.LATTICE # list of section keys in a Lattice z.B. [lE,HE,...], in order from left=entance to right=exut
        DEBUG_OFF(sectionIDs)
        DEBUG_OFF('make_lattice for sollteilchen\n'+util.PARAMS['sollteilchen'].string())

        for sectionID in sectionIDs:
            elementIDs = sections.get(sectionID)
            for elementID in elementIDs:
                element        = util.ELEMENTS.get(elementID)
                """add sectionID and elementID"""
                element['sec'] = sectionID 
                element['ID']  = elementID 
                item           = {elementID:element}  # repack {ID:{attributes}} for instanciate_element(...)
                """INSTANCIATE ELM._Node objects"""
                instance = instanciate_element(item)
                if isinstance(instance,ELM._Node):
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

    results = namedtuple('InputParseResult',['SECTIONSx','SECTIONSu','LATTICE','FLAGS','PARAMETERS','ELEMENTS'])
    results = parse(in_data)

    flags = proces_flags(results.FLAGS)
    DEBUG_OFF('global FLAGS after proces_flags():')
    DEBUG_OFF(util.FLAGS)

    parameters = proces_parameters(results.PARAMETERS)
    DEBUG_OFF('global PARAMS after proces_parameters():')
    DEBUG_OFF(util.PARAMS)
    
    elements = proces_elements(results.ELEMENTS)
    DEBUG_OFF('ELEMENTS after proces_elements():')
    DEBUG_OFF(util.ELEMENTS)

    sections = proces_sections(results)
    DEBUG_OFF('SECTIONS after proces_sections():')
    DEBUG_OFF(util.SECTIONS)

    sectionIDs = proces_lattice(results.LATTICE)
    DEBUG_OFF('LATTICE after proces_lattice():')
    DEBUG_OFF(util.LATTICE)

    util.PARAMS['sollteilchen'](tkin=util.PARAMS['injection_energy'])

    lattice = make_lattice()
    DEBUG_OFF('lattice_generator >>{}'.format(lattice.string()))
    util.SUMMARY['lattice length [m]'] = util.PARAMS['lattice_length']  = lattice.length
    DEBUG_OFF('SUMMARY in factory() {}'.format(util.SUMMARY))
    # end of factory(...)
    return lattice

def test0(input_file):
    print('---------------------------------TEST0')
    wfl= []
    fileobject=open(input_file,'r')
    wfl= yaml.load(fileobject)
    print(yaml.dump(wfl,default_flow_style=True))
    for i,v in iter(wfl.items()):
        print(i,' =\t',v)
    seg = wfl['SEGMENTS']
    print(seg)
    print('=== segment ===')
    for i in seg:
        print(i)
    lattice = wfl['LATTICE']
    print('=== lattice ===')
    for l in lattice:
        # for i in l:
        #     print(i)
        print(l)
def test1(input_file):
    print('---------------------------------TEST1')
    lattice = factory_new(input_file)
    print('%%%%%%%%%%%%%%%%%%% Left------->Right')
    for cnt,node in enumerate(iter(lattice)):
        print(cnt,'{:38s} {:38s}'.format(repr(node),repr(node.next)))
    print('\n%%%%%%%%%%%%%%%%%%% Right------->Left')
    lattice.toggle_iteration()
    for cnt,node in enumerate(iter(lattice)):
        print(cnt,'{:38s} {:38s}'.format(repr(node),repr(node.prev)))

if __name__ == '__main__':
    input_file = 'yml/tmpl_25.10.2021_new.yml'
    # input_file = 'yml/simuIN.yml'
    # input_file = 'nwlat/nwlatIN.yml'
    # test0(input_file)
    test1(input_file)

