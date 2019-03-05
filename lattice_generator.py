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

from setutil import PARAMS,FLAGS,SUMMARY,DEBUG,DEBUG_ON,DEBUG_OFF,Twiss
import elements as ELM
from lattice import Lattice
from Ez0 import SFdata,V0
import marker_actions as MRK

DEBUG_MODULE = DEBUG_OFF

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
def flatten(lis):
    """Given a list, possibly nested to any level, return it flattened."""
    new_lis = []
    for item in lis:
        if type(item) == type([]):
            new_lis.extend(flatten(item))
        else:
            new_lis.append(item)
    return new_lis

def liofd2d(l):
    """Transform list of dicts to dict"""
    return {k:v for d in l for k,v in d.items()}

def replace_QF_with_QFth_lattice(slices,k0,length,label,particle,aperture):
    lattice = Lattice()
    thinlen = length/slices
    for nb in range(slices):
        if FLAGS['express']:
            instance = ELM.QFthx(k0=k0,length=thinlen,label=label,particle=particle,aperture=aperture)
        else:
            instance = ELM.QFth(k0=k0,length=thinlen,label=label,particle=particle,aperture=aperture)
        lattice.add_element(instance)
    return lattice

def replace_QD_with_QDth_lattice(slices,k0,length,label,particle,aperture):
    lattice = Lattice()
    thinlen = length/slices
    for nb in range(slices):
        if FLAGS['express']:
            instance = ELM.QDthx(k0=k0,length=thinlen,label=label,particle=particle,aperture=aperture)
        else:
            instance = ELM.QDth(k0=k0,length=thinlen,label=label,particle=particle,aperture=aperture)
        lattice.add_element(instance)
    return lattice
def instanciate_element(item):
    DEBUG_MODULE('instanciate_element: instanciate {}'.format(item))
    key = item[0]
    attributes = item[1]
    # aperture   = PARAMS['aperture']    # default aperture
    if key == 'D':
        label    = attributes['ID']
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
        kq       = dBdz/PARAMS['sollteilchen'].brho
        Bpole    = dBdz*aperture
        if slices > 1:
            instance = replace_QF_with_QFth_lattice(slices,kq,length,label,PARAMS['sollteilchen'],aperture)
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
        kq       = dBdz/PARAMS['sollteilchen'].brho
        Bpole    = dBdz*aperture
        if slices > 1:
            instance = replace_QD_with_QDth_lattice(slices,kq,length,label,PARAMS['sollteilchen'],aperture)
        elif slices <= 1:
            instance = ELM.QD(k0=kq,length=length,label=label,aperture=aperture)
        instance['label']  = label
        instance['dBdz']   = dBdz
        instance['bore']   = aperture
        instance['Bpole']  = Bpole

    elif key == 'RFG':
        label     = attributes['ID']
        PhiSoll   = radians(get_mandatory(attributes,"PhiSync",label))
        fRF       = get_mandatory(attributes,"fRF",label)
        gap       = get_mandatory(attributes,'gap',label)
        aperture  = get_mandatory(attributes,'aperture',label)
        dWf       = FLAGS['dWf']
        mapping   = PARAMS['mapping']
        if mapping == None:
            mapping = 't3d'
        if mapping == 'ttf' or mapping == 'dyn' or mapping == 'oxal': # SF-data
            fname     = get_mandatory(attributes,"SFdata",label)
            EzPeak    = get_mandatory(attributes,"EzPeak",label)
            if fname not in PARAMS:
                PARAMS[fname] = SFdata(fname,EzPeak=EzPeak)
            EzAvg = PARAMS[fname].EzAvg
            instance  =  ELM.RFG(EzAvg=EzAvg,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,mapping=mapping,dWf=dWf,aperture=aperture,SFdata=PARAMS[fname])
            instance['EzPeak'] = EzPeak
            pass
        else:
            EzAvg     = get_mandatory(attributes,"EzAvg",label)
            instance  = ELM.RFG(EzAvg=EzAvg,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,mapping=mapping,dWf=dWf,aperture=aperture)
            instance['EzPeak'] = None
        instance['label']    = label
        instance['PhiSoll']  = PhiSoll
        instance['fRF']      = fRF
        instance['gap']      = gap
        instance['aperture'] = aperture
        instance['dWf']      = dWf
        instance['mapping']  = mapping

    elif key == 'RFC':
        label     = attributes['ID']
        PhiSoll   = radians(get_mandatory(attributes,"PhiSync",label))
        fRF       = get_mandatory(attributes,"fRF",label)
        gap       = get_mandatory(attributes,'gap',label)
        aperture  = get_mandatory(attributes,'aperture',label)
        dWf       = FLAGS['dWf']
        length    = get_mandatory(attributes,'length',label)
        mapping   = PARAMS['mapping']
        if mapping == None:
            mapping = 't3d'
        if mapping == 'ttf' or mapping == 'dyn' or mapping == 'oxal': # SF-data
            fname     = get_mandatory(attributes,"SFdata",label)
            EzPeak    = get_mandatory(attributes,"EzPeak",label)
            if fname not in PARAMS:
                PARAMS[fname] = SFdata(fname,EzPeak=EzPeak)
            EzAvg = PARAMS[fname].EzAvg
            instance  =  ELM.RFC(EzAvg=EzAvg,label=label,PhiSoll=PhiSoll,fRF=fRF,gap=gap,aperture=aperture,dWf=dWf,length=length,mapping=mapping,SFdata=PARAMS[fname])
            instance['EzPeak'] = EzPeak
            pass
        else:
            EzAvg     = get_mandatory(attributes,"EzAvg",label)
            instance  =  ELM.RFC(EzAvg=EzAvg,label=label,PhiSoll=PhiSoll,fRF=fRF,gap=gap,aperture=aperture,dWf=dWf,length=length,mapping=mapping)
            instance['EzPeak'] = None
        instance['label']    = label
        instance['PhiSoll']  = PhiSoll
        instance['fRF']      = fRF
        instance['gap']      = gap
        instance['aperture'] = aperture
        instance['dWf']      = dWf
        instance['length']   = length
        instance['mapping']  = mapping

    elif key == 'GAP':
        label     = attributes['ID']
        gap       = get_mandatory(attributes,'gap',label)
        EzAvg     = get_mandatory(attributes,"EzAvg",label)
        PhiSoll   = radians(get_mandatory(attributes,"PhiSync",label))
        fRF       = get_mandatory(attributes,"fRF",label)
        dWf       = FLAGS['dWf']
        aperture  = get_mandatory(attributes,'aperture',label)
        instance  =  ELM.GAP(EzAvg=EzAvg,PhiSoll=PhiSoll,fRF=fRF,label=label,gap=gap,dWf=dWf,aperture=aperture)
        instance['label']   = label
        instance['gap']     = gap
        instance['EzAvg']   = EzAvg
        instance['PhiSoll'] = PhiSoll
        instance['fRF']     = fRF
        instance['dWf']     = dWf

    elif key == 'MRK':
        label     = attributes['ID']
        action    = get_mandatory(attributes,'action',label)
        if 'scatter' == action:
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
    return (label,instance)

def factory(input_file):
    """ factory creates a lattice from input-file """
#--------
    def read_elements(in_data):
        element_list = in_data['ELEMENTS']
        elements = liofd2d(element_list)
        # add {'ID':key} to attribute list
        IDs = elements.keys()
        for ID in IDs:
            attList = elements[ID]
            attList.append({'ID':ID})
        return elements
# --------
    def read_flags(in_data):
        """returns a dict of flags"""
        flags_list = in_data['FLAGS']
        flags = liofd2d(flags_list) if flags_list != None else {}
        if 'accON' in flags:
            if flags['accON']:
                FLAGS['dWf'] = 1.
                SUMMARY['accON'] = True
            else:
                FLAGS['dWf'] = 0.
                SUMMARY['accON'] = False
        if 'periodic'    in flags: FLAGS['periodic'] = flags['periodic']
        if 'egf'         in flags: FLAGS['egf']      = flags['egf']
        if 'sigma'       in flags: FLAGS['sigma']    = flags['sigma']
        if 'KVout'       in flags: FLAGS['KVout']    = flags['KVout']
        if 'verbose'     in flags: FLAGS['verbose']  = flags['verbose']
        if 'express'     in flags: FLAGS['express']  = flags['express']
        if 'useaper'     in flags: FLAGS['useaper']  = flags['useaper']
        if 'bucket'      in flags: FLAGS['bucket']   = flags['bucket']
        if 'csTrak'      in flags: FLAGS['csTrak']   = flags['csTrak']
        if 'pspace'      in flags: FLAGS['pspace']   = flags['pspace']
        return flags
# --------
    def read_sections(in_data):
        """ returns a list of section names """
        sec_list = []
        use_sections = True
        try:     ## can fail because sections are not mandatory
            sec_list = in_data['SECTIONS']
            sec_list = flatten(sec_list)
        except:
            use_sections = False
        PARAMS['sections'] = sec_list
        FLAGS['sections']  = use_sections
        return sec_list
# --------
    def read_parameters(in_data):
        """ returns a dict of parameters """
        parameter_list = in_data['PARAMETERS']
        parameters     = liofd2d(parameter_list)
        if 'frequency'        in parameters: PARAMS['frequenz']         = parameters['frequency']
        if 'B_grad_f'         in parameters: PARAMS['qf_gradient']      = parameters['B_grad_f']
        if 'B_grad_d'         in parameters: PARAMS['qd_gradient']      = parameters['B_grad_d']
        if 'Tkin'             in parameters: PARAMS['injection_energy'] = parameters['Tkin']
        if 'EzAvg'            in parameters: PARAMS['EzAvg']            = parameters['EzAvg']
        if 'phi_sync'         in parameters: PARAMS['phisoll']          = parameters['phi_sync']
        if 'gap'              in parameters: PARAMS['gap']              = parameters['gap']
        if 'cav_len'          in parameters: PARAMS['cavity_laenge']    = parameters['cav_len']
        if 'ql'               in parameters: PARAMS['ql']               = parameters['ql']
        if 'windings'         in parameters: PARAMS['nbwindgs']         = parameters['windings']
        if 'nbsigma'          in parameters: PARAMS['nbsigma']          = parameters['nbsigma']
        if 'aperture'         in parameters: PARAMS['aperture']         = parameters['aperture'] 
        if 'emitx_i'          in parameters: PARAMS['emitx_i']          = parameters['emitx_i']
        if 'emity_i'          in parameters: PARAMS['emity_i']          = parameters['emity_i']
        if 'betax_i'          in parameters: PARAMS['betax_i']          = parameters['betax_i']
        if 'betay_i'          in parameters: PARAMS['betay_i']          = parameters['betay_i']
        if 'alfax_i'          in parameters: PARAMS['alfax_i']          = parameters['alfax_i']
        if 'alfay_i'          in parameters: PARAMS['alfay_i']          = parameters['alfay_i']
        if 'mapping'          in parameters: PARAMS['mapping']          = parameters['mapping']
        if 'DT2T'             in parameters: PARAMS['DT2T']             = parameters['DT2T']
        if 'lattvers'         in parameters: PARAMS['lattice_version']  = parameters['lattvers']

        PARAMS['lamb']           = PARAMS['lichtgeschwindigkeit']/PARAMS['frequenz']
        PARAMS['wellenlänge']    = PARAMS['lichtgeschwindigkeit']/PARAMS['frequenz']
        PARAMS['spalt_spannung'] = PARAMS['EzAvg']*PARAMS['gap']
        return parameters
#--------
    def get_flattened_lattice_list(in_data):
        """ read_and_flatten lattice from (in_data)."""
        lattice_list = in_data['LATTICE']
        lattice_list = flatten(lattice_list)
        N = lattice_list[0]   # Duplikator
        plist = lattice_list[1:]
        qlist = plist.copy()
        for i in range(N-1):
            qlist += plist
        lattice_list=qlist
        return lattice_list
#--------
    def make_lattice(latticeList,in_data):
        """ instanciate all elements from flattened node list"""
        lattice = Lattice()
        DEBUG_OFF('make_lattice for sollteilchen\n'+PARAMS['sollteilchen'].string())
        elements = read_elements(in_data)
        for ID in lattice_list:
            element      = elements[ID]
            element      = liofd2d(element)
            elementClass = element['type']
            elmItem      = (elementClass,element)
            # !!INSTANCIATE!!
            (label,instance) = instanciate_element(elmItem)
            section = instance.section if FLAGS['sections'] else '*'
            DEBUG_MODULE('instance {} {} {}'.format(label,instance,section))
            # add element instance to lattice
            if isinstance(instance,ELM._Node):
                lattice.add_element(instance)
        #     elif isinstance(instance,Lattice):
        #         lattice.concat(instance)       # concatenate partial with lattice
        return lattice   # the complete lattice

## factory body --------
    SUMMARY['input file'] = PARAMS['input_file'] = input_file

    with open(input_file,'r') as fileobject:
        try:
            in_data = yaml.load(fileobject)
        except Exception as ex:
            warnings.showwarning(
                    'InputError: {} - STOP'.format(str(ex)),
                    UserWarning,
                    'lattice_generator.py',
                    'factory()',
                    )
            sys.exit(1)
    fileobject.close()

    read_flags(in_data)
    read_sections(in_data)
    read_parameters(in_data)
    DEBUG_MODULE('PARAMS after read_parameters()',PARAMS)
    
    # lattice_list is a flat list of node IDs
    lattice_list = get_flattened_lattice_list(in_data)
    DEBUG_MODULE('latticeList in factory()',lattice_list)
    # __call__ sollteilchen energy
    PARAMS['sollteilchen'](tkin=PARAMS['injection_energy'])
    lattice = make_lattice(lattice_list,in_data)
    # DEBUG_MODULE('lattice_generator >>\n',lattice.string())
    SUMMARY['lattice length [m]'] = PARAMS['lattice_length']  = lattice.length
    DEBUG_OFF('SUMMARY in factory()',SUMMARY)
    
    # end of factory(...)
    return lattice

def parse_and_fabric(input_file):   # delegates to factory
    return factory(input_file)

# utilities
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
        for i in l:
            print(i)

def test1(input_file):
    print('---------------------------------TEST1')
    lattice = parse_and_fabric(input_file)
    print('%%%%%%%%%%%%%%%%%%% Left------->Right')
    for cnt,node in enumerate(iter(lattice)):
        print(cnt,'{:38s} {:38s}'.format(repr(node),repr(node.next)))
    print('\n%%%%%%%%%%%%%%%%%%% Right------->Left')
    lattice.toggle_iteration()
    for cnt,node in enumerate(iter(lattice)):
        print(cnt,'{:38s} {:38s}'.format(repr(node),repr(node.prev)))

def test2(input_file):
    print('---------------------------------TEST2')
    with open(input_file,'r') as fileobject:
        try:
            indat = yaml.load(fileobject)
        except Exception as ex:
            warnings.showwarning(
                    'InputError: {} - STOP'.format(str(ex)),
                    UserWarning,
                    'lattice_generator.py',
                    'factory()',
                    )
            sys.exit(1)
    fileobject.close()
    ky = indat.keys()
    for k in ky:
        print()
        print(k)
        klist = indat[k]
        print(klist)
        if not klist: continue
        nlist = flatten(klist)
        if k == 'LATTICE':
            N = nlist[0]
            plist = nlist[1:]
            qlist = plist.copy()
            for i in range(N-1):
                qlist += plist
            nlist=qlist
        print(nlist)

if __name__ == '__main__':
    input_file = 'yml/simuINstat.yml'
    # input_file = 'nwlat/nwlatIN.yml'
    test0(input_file)
    # test1(input_file)
    test2(input_file)

