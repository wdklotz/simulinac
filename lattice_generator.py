#!/Users/klotz/SIMULINAC_env/bin/python
# -*- coding: utf-8 -*-
__version__='11.0.2.4'
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
import sys,os
import yaml
import warnings
import unittest
import setutil             as UTIL
import elements            as ELM
import PsMarkerAgent       as PSMKR
import PoincareMarkerAgent as PCMKR
import lattice_parser_2    as LP2
import lattice             as LAT
import Ez0                 as EZ
from math    import radians
from T3D_M   import T3D_G
from OXAL_M  import OXAL_G
from Base_M  import Base_G
from TTF_M   import TTF_G

wrapRED = UTIL.wrapRED

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
        raise(UserWarning(wrapRED('Mandatory attribute "{}" missing for element "{}"'.format(key,item))))
        sys.exit(1)
    return res

def mandatory_warning(key,item):
    raise(UserWarning(wrapRED('Mandatory attribute "{}" missing for element "{}"'.format(key,item))))
    sys.exit(1)

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
        elif type == 'SD':
            alpha          = get_mandatory(attributes,'alpha',ID)
            rho            = get_mandatory(attributes,'rho',ID)
            aperture       = attributes.get('aperture',None)
            instance       = ELM.SD(ID, alpha, rho, aperture=aperture)
            ELEMENT['sec'] = attributes.get('sec','?')
            ELEMENT['B0 [T]'] = 3.3356/rho*UTIL.PARAMS['injection_energy']*1e-3
        elif type == 'RD':
            alpha          = get_mandatory(attributes,'alpha',ID)
            rho            = get_mandatory(attributes,'rho',ID)
            aperture       = attributes.get('aperture',None)
            t3d_wedge      = attributes.get('t3d_wedge',False)
            wedge          = ELM.Wedge(alpha/2., rho, t3d_wedge=t3d_wedge)
            instance       = ELM.RD(ID, alpha, rho, wedge, aperture=aperture)
            ELEMENT['sec'] = attributes.get('sec','?')
            ELEMENT['length'] = instance.length
            ELEMENT['B0 [T]'] = 3.3356/rho*UTIL.PARAMS['injection_energy']*1e-3
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
        elif type == 'RFG_OLD':
            phiSoll          = M.radians(get_mandatory(attributes,"PhiSync",ID))
            freq             = float(get_mandatory(attributes,"freq",ID))
            aperture         = get_mandatory(attributes,'aperture',ID)
            particle         = UTIL.Proton(UTIL.PARAMS['injection_energy'])
            position         = (0.,0.,0.)
            dWf              = UTIL.FLAGS['dWf']
            ELEMENT['sec']   = attributes.get('sec','?')

            mapping = UTIL.FLAGS.get('mapping')    # global FLAG overrides local mapping
            if mapping == None:
                mapping = attributes.get('mapping','t3d')
            ELEMENT['mapping'] = mapping   # maybe overriden by global mapping

            if mapping == 't3d'or mapping == "simple":
                fieldtab = attributes.get('SFdata',None)
                if fieldtab == None:
                    EzPeak    = get_mandatory(attributes,"EzPeak",ID)
                    gap       = get_mandatory(attributes,'gap',ID)
                    cavlen    = get_mandatory(attributes,'cavlen',ID)
                    SFdata    = None
                    instance  = ELM.RFG(ID,EzPeak=EzPeak,phisoll=phiSoll,gap=gap,cavlen=cavlen,freq=freq,SFdata=SFdata,particle=particle,position=position,aperture=aperture,dWf=dWf,mapping=mapping)
                    instance.gap_object = instance   # instance gets reference to itself as member
                    ELEMENT['HE_Gap'] ='ignored'
                    pass
                else:
                    EzPeak    = get_mandatory(attributes,"EzPeak",ID)
                    cavlen    = get_mandatory(attributes,'cavlen',ID)  # cavity length [m]
                    HE_Gap    = get_mandatory(attributes,'HE_Gap',ID)  # hard edge gap [m]
                    sfdata    = EZ.SFdata.InstanciateAndScale(fieldtab,EzPeak=EzPeak,L=cavlen/2.*100.)   # scaled field distribution
                    (dummy,HE_EzPeak) = sfdata.hardEdge(HE_Gap*100.)   # hard edge: HE_Gap [cm], HE_EzPeak [MV/m]
                    instance  = ELM.RFG(ID,EzPeak=HE_EzPeak,phisoll=phiSoll,gap=HE_Gap,cavlen=cavlen,freq=freq,SFdata=sfdata,particle=particle,position=position,aperture=aperture,dWf=dWf,mapping=mapping)
                    instance.gap_object = instance   # instance gets reference to itself as member
                    ELEMENT['HE_EzPeak'] = HE_EzPeak
                    ELEMENT['gap']       = 'ignored'
                    pass
            elif mapping == 'oxal':
                fieldtab = get_mandatory(attributes,"SFdata",ID)
                EzPeak   = get_mandatory(attributes,"EzPeak",ID)
                cavlen   = get_mandatory(attributes,"cavlen",ID)
                SFdata   = EZ.SFdata.InstanciateAndScale(fieldtab,EzPeak=EzPeak,L=cavlen/2.*100.)   # scaled field distribution
                oxa      = OXA.OXAL_G(ID,EzPeak=EzPeak,phisoll=phiSoll,cavlen=cavlen,freq=freq,SFdata=SFdata,particle=particle,position=position,aperture=aperture,dWf=dWf)
                instance = ELM.RFG(ID,particle=particle,freq=freq,mapping=mapping)
                instance.gap_object = oxa   # instance gets reference to OXAL_G as member
                ELEMENT['HE_Gap'] ='ignored'
                ELEMENT['gap']    ='ignored'
                pass
            elif mapping == 'ttf':
                raise(UserWarning(wrapRED(f'Mapping "{mapping}" is not available')))
                sys.exit()
            else:
                raise(UserWarning(wrapRED(f'Mapping "{mapping}" is not available')))
                sys.exit()

            # instance created with parameters, configure with matrix
            instance.dispatch_model_matrix()
            pass
        elif type == 'RFG':
            aperture   = attributes.get('aperture')
            cavlen     = attributes.get('cavlen')
            EzPeak     = attributes.get('EzPeak')
            freq       = attributes.get('freq')
            gap        = attributes.get('gap')
            HE_Gap     = attributes.get('HE_Gap')
            mapping    = ELEMENT['mapping'] = UTIL.FLAGS.get('mapping')
            phisoll    = attributes.get('phisoll')
            SFdata     = attributes.get('SFdata')
            sec        = attributes.get('sec','?')
            
            if mapping   == 't3d':
                if SFdata == None:
                    if(aperture == None): mandatory_warning("aperture",ID) # [MV/m] requested aperture
                    if(EzPeak   == None): mandatory_warning("EzPeak",ID)   # [MV/m] requested peak field
                    if(freq     == None): mandatory_warning("freq",ID)   # [MV/m] requested frequency
                    if(gap      == None): mandatory_warning("gap",ID)      # [MV/m] requested gap
                    if(phisoll  == None): mandatory_warning("phisoll",ID)  # [MV/m] requested synch phase

                    gap_parameters = dict(
                        aperture  = aperture,          # [m]
                        cavlen    = cavlen,            # [m]
                        EzPeak    = EzPeak,            # [MV/m]
                        freq      = freq,              # [Hz] 
                        gap       = gap,               # [m]
                        HE_Gap    = HE_Gap,               # [m]
                        phisoll   = radians(phisoll),  # [radians] 
                        sec       = sec,               # string
                        SFdata    = SFdata,
                    )
                    instance = ELM.RFG(ID)
                    instance.register(T3D_G())         # register T3D
                    instance.configure(**gap_parameters)
                    if cavlen != None: ELEMENT['cavlen'] ='ignored'
                    if HE_Gap != None: ELEMENT['HE_Gap'] ='ignored'
                    if SFdata != None: ELEMENT['SFdata'] ='ignored'
                    pass
                elif SFdata != None:
                    if(aperture == None): mandatory_warning("aperture",ID)
                    if(cavlen   == None): mandatory_warning("cavlen",ID)
                    if(EzPeak   == None): mandatory_warning("EzPeak",ID)
                    if(freq     == None): mandatory_warning("EzPeak",ID)
                    if(HE_Gap   == None): mandatory_warning("HE_Gap",ID)
                    if(phisoll  == None): mandatory_warning("phisoll",ID)

                    # fieldtab = attributes.get('SFdata',None)
                    sfdata = EZ.SFdata.InstanciateAndScale(SFdata,EzPeak=EzPeak,L=cavlen*100.)   # scaled field distribution
                    (dummy,HE_EzPeak) = sfdata.hardEdge(HE_Gap*100)     # [MV/m] equivalent hard edge peak field
                    ELEMENT['HE_EzPeak'] = HE_EzPeak
                    gap_parameters = dict(
                        aperture  = aperture,
                        cavlen    = cavlen,
                        EzPeak    = HE_EzPeak,
                        freq      = freq,       # [Hz]  requested RF frequenz
                        gap       = HE_Gap,
                        HE_Gap    = HE_Gap,               # [m]
                        phisoll   = radians(phisoll),    # [radians] requested soll phase
                        sec       = sec,
                        SFdata    = sfdata,
                    )
                    instance = ELM.RFG(ID)
                    instance.register(T3D_G())
                    instance.configure(**gap_parameters)
                    if gap != None: ELEMENT['gap'] = 'ignored'
                    pass
            elif mapping == 'simple':
                raise(UserWarning(wrapRED(f'Mapping "{mapping}" is not ready. Must be implemented')))
                sys.exit()
            elif mapping == 'oxal':
                if(aperture == None): mandatory_warning("aperture",ID) # [MV/m] requested aperture
                if(cavlen   == None): mandatory_warning("cavlen",ID)
                if(EzPeak   == None): mandatory_warning("EzPeak",ID)
                if(freq     == None): mandatory_warning("freq",ID)   # [MV/m] requested frequency
                if(phisoll  == None): mandatory_warning("phisoll",ID)  # [MV/m] requested synch phase
                if(SFdata   == None): mandatory_warning("SFdata",ID)

                sfdata = EZ.SFdata.InstanciateAndScale(SFdata,EzPeak=EzPeak,L=cavlen*100.)   # scaled field distribution
                gap_parameters = dict(
                    aperture  = aperture,
                    cavlen    = cavlen,
                    EzPeak    = EzPeak,
                    freq      = freq,
                    gap       = gap,
                    HE_Gap    = HE_Gap,
                    phisoll   = radians(phisoll),
                    sec       = sec,
                    SFdata    = sfdata,
                )
                instance = ELM.RFG(ID)
                instance.register(OXAL_G())
                instance.configure(**gap_parameters)
                if HE_Gap != None: ELEMENT['HE_Gap'] ='ignored'
                if gap    != None: ELEMENT['gap']    ='ignored'
                pass
            elif mapping == 'base':
                if(aperture == None): mandatory_warning("aperture",ID) # [MV/m] requested aperture
                if(EzPeak   == None): mandatory_warning("EzPeak",ID)
                if(freq     == None): mandatory_warning("freq",ID)   # [MV/m] requested frequency
                if(gap      == None): mandatory_warning("gap",ID)
                if(phisoll  == None): mandatory_warning("phisoll",ID)  # [MV/m] requested synch phase
                gap_parameters = dict(
                    aperture  = aperture,
                    cavlen    = cavlen,
                    EzPeak    = EzPeak,
                    freq      = freq,
                    gap       = gap,
                    HE_Gap    = HE_Gap,
                    phisoll   = radians(phisoll),
                    sec       = sec,
                    SFdata    = SFdata,
                )
                instance = ELM.RFG(ID)
                instance.register(Base_G())
                instance.configure(**gap_parameters)
                if HE_Gap != None: ELEMENT['HE_Gap'] = 'ignored'
                if cavlen != None: ELEMENT['cavlen'] = 'ignored'
                if SFdata != None: ELEMENT['SFdata'] = 'ignored'
                pass
            elif mapping == 'ttf':
                if(aperture == None): mandatory_warning("aperture",ID) # [MV/m] requested aperture
                if(cavlen == None): mandatory_warning("cavlen",ID) # [MV/m] requested aperture
                if(EzPeak   == None): mandatory_warning("EzPeak",ID)
                if(freq     == None): mandatory_warning("freq",ID)   # [MV/m] requested frequency
                if(phisoll  == None): mandatory_warning("phisoll",ID)  # [MV/m] requested synch phase
                if(SFdata   == None): mandatory_warning("SFdata",ID)  # [MV/m] requested synch phase

                sfdata = EZ.SFdata.InstanciateAndScale(SFdata,EzPeak=EzPeak,L=cavlen/2.*100.)   # scaled field distribution
                gap_parameters = dict(
                    aperture  = aperture,
                    cavlen    = cavlen,
                    EzPeak    = EzPeak,
                    freq      = freq,
                    gap       = gap,
                    HE_Gap    = HE_Gap,
                    phisoll   = radians(phisoll),
                    sec       = sec,
                    SFdata    = sfdata,
                )
                instance = ELM.RFG(ID)
                instance.register(TTF_G())
                instance.configure(**gap_parameters)
                if gap    != None: ELEMENT['gap']    ='ignored'
                if HE_Gap != None: ELEMENT['HE_Gap'] ='ignored'
                pass
            elif mapping == 'dyn':
                raise(UserWarning(wrapRED(f'Mapping "{mapping}" is not available any more')))
                sys.exit()
            else:
                raise(UserWarning(wrapRED(f'Invalid mapping "{mapping}"')))
                sys.exit()
        elif type == 'RFC':
            raise(UserWarning(wrapRED(f'Element "{type}" is not ready. Must be implemented')))
            sys.exit()
            phisoll          = M.radians(get_mandatory(attributes,"PhiSync",ID))
            freq             = float(get_mandatory(attributes,"freq",ID))
            aperture         = get_mandatory(attributes,'aperture',ID)
            particle         = UTIL.Proton(UTIL.PARAMS['injection_energy'])
            position         = (0.,0.,0.)
            dWf              = UTIL.FLAGS['dWf']
            ELEMENT['sec']   = attributes.get('sec','?')

            mapping = UTIL.FLAGS.get('mapping')    # global mapping FLAG overrides individual mapping
            if mapping == None:
                mapping = attributes.get('mapping','t3d')
            ELEMENT['mapping'] = mapping  
            # if mapping == 't3d' or mapping == "simple":
            if mapping == 't3d':
                fieldtab = attributes.get('SFdata',None)
                if fieldtab == None:
                    EzPeak    = get_mandatory(attributes,"EzPeak",ID)   # [MV/m] requested peak field
                    gap       = get_mandatory(attributes,'gap',ID)      # [m] requested gap length
                    cavlen    = get_mandatory(attributes,'cavlen',ID)   # [m] requested cavity length
                    gap_parameters = dict(
                        EzPeak    = EzPeak,
                        phisoll   = phisoll,         # [radians] requested soll phase
                        gap       = gap,
                        cavlen    = cavlen,
                        freq      = freq,            # [Hz]  requested RF frequenz
                        particle  = particle,
                        position  = position,
                        aperture  = aperture
                    )
                    # instance  = ELM.RFC(ID,EzPeak,phiSoll,gap,cavlen,freq,SFdata=0,particle=particle,position=position,aperture=aperture,dWf=dWf,mapping=mapping)
                    instance = ELM.RFG(ID)
                    instance.register_mapping(T3D_G())
                    instance.configure(**gap_parameters)
                    ELEMENT['HE_Gap'] ='ignored'
                    pass
                else:
                    EzPeak    = get_mandatory(attributes,"EzPeak",ID)   # [MV/m] requested peak field
                    cavlen    = get_mandatory(attributes,'cavlen',ID)   # [m] requested cavity length
                    HE_Gap    = get_mandatory(attributes,'HE_Gap',ID)   # [m] requested gap length [m]
                    sfdata    = EZ.SFdata.InstanciateAndScale(fieldtab,EzPeak=EzPeak,L=cavlen/2*100.)   # scaled field distribution
                    (dummy,HE_EzPeak) = sfdata.hardEdge(HE_Gap*100)     # [MV/m] equivalent hard edge peak field
                    gap_parameters = dict(
                        EzPeak    = HE_EzPeak,
                        phisoll   = phisoll,    # [radians] requested soll phase
                        gap       = HE_Gap,
                        cavlen    = cavlen,
                        freq      = freq,       # [Hz]  requested RF frequenz
                        particle  = particle,
                        position  = position,
                        aperture  = aperture
                    )
                    # instance  = ELM.RFC(ID,HE_EzPeak,phiSoll,HE_Gap,cavlen,freq,SFdata=sfdata,particle=particle,position=position,aperture=aperture,dWf=dWf,mapping=mapping)
                    instance = ELM.RFG(ID)
                    instance.register_mapping(T3D_G())
                    instance.configure(**gap_parameters)
                    ELEMENT['HE_EzPeak'] = HE_EzPeak
                    ELEMENT['gap']       = 'ignored'
                    pass
            elif mapping == 'oxal':
                fieldtab  = get_mandatory(attributes,'SFdata',ID)
                EzPeak    = get_mandatory(attributes,"EzPeak",ID)
                cavlen    = get_mandatory(attributes,'cavlen',ID)  # cavity length [m]
                sfdata    = EZ.SFdata.InstanciateAndScale(fieldtab,EzPeak=EzPeak,L=cavlen/2.*100.)   # scaled field distribution
                instance  = ELM.RFC(ID,EzPeak,phiSoll,0.,cavlen,freq,SFdata=sfdata,particle=particle,position=position,aperture=aperture,dWf=dWf,mapping=mapping)
                ELEMENT['HE_Gap'] ='ignored'
                ELEMENT['gap']    ='ignored'
                pass
            elif mapping == 'ttf':
                raise(UserWarning(wrapRED(f'Mapping "{mapping}" is not ready. Must be implemented')))
                sys.exit()
            elif mapping == 'dyn':
                raise(UserWarning(wrapRED(f'Mapping "{mapping}" is not available any more')))
                sys.exit()
            else:
                raise(UserWarning(wrapRED(f'Invalid mapping "{mapping}"')))
                sys.exit()
        elif type == 'GAP':
            gap       = get_mandatory(attributes,'gap',ID)
            EzPeak    = get_mandatory(attributes,"EzPeak",ID)
            phiSoll   = M.radians(get_mandatory(attributes,"PhiSync",ID))
            freq      = float(get_mandatory(attributes,"freq",ID))
            dWf       = UTIL.FLAGS['dWf']
            aperture  = get_mandatory(attributes,'aperture',ID)
            EzAvg     = EzPeak
            instance  =  ELM.GAP(ID,EzAvg,phiSoll,gap,freq,aperture=aperture,dWf=dWf)
            ELEMENT['EzPeak'] = EzPeak
            ELEMENT['sec']    = attributes.get('sec','?')
            ELEMENT['EzAvg']  = EzAvg
        elif type == 'MRK':
            # active  = attributes.get('active',UTIL.FLAGS['maction'])
            active  = get_mandatory(attributes,'active',ID)
            action  = get_mandatory(attributes,'action',ID)
            viseo   = attributes.get('viseo',3)
            ELEMENT = ELEMENT
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
                raise(UserWarning(wrapRED('Unknown marker ACTION encountered: "{}"'.format(action))))
                sys.exit(1)
        else:
            raise(UserWarning(wrapRED('Unknown element "{}" encountered.'.format(type))))
            sys.exit(1)
    return instance

def factory(input_file,stop=None):
    """ factory reads FLAGS, PARAMETERS and creates a lattice from yaml input """
    
    def process_flags(flags):
        """ external FLAGs """        
        res = dict(
            accON    = flags.get('accON',True),               # acceleration ON
            periodic = flags.get('periodic',False),           # periodic lattice? default
            egf      = flags.get('egf',False),                # emittance grow flag default
            sigma    = flags.get('sigma',True),               # beam sizes by sigma-tracking
            KVout    = flags.get('KVout',False),              # print a dictionary of Key-Value pairs, no display
            GDisp    = flags.get('GDisp',True),               # display, show graphics, default show
            verbose  = flags.get('verbose',0),                # print flag default = 0
            useaper  = flags.get('useaper',False),            # use aperture check for quads and rf-gaps
            bucket   = flags.get('bucket',False),             # plot bucket
            csTrak   = flags.get('csTrak',True),              # plot CS trajectories
            maction  = flags.get('maction',True),             # call marker actions
            envelope = flags.get('envelope',False),           # plot transverse envelopes
            mapping  = flags.get('mapping','t3d')             # global mapping overrides individula mapping (default to t3d)
        )
        """ internal FLAGs """        
        res['dWf'] = 1 if res['accON'] else 0    # acceleration on/off flag 1=on,0=off
        res['non_linear_mapping'] = False
        if res['mapping'] == 'ttf' or res['mapping'] == 'base' or res['mapping'] == 'dyn': 
            res['non_linear_mapping'] = True
        return res

    def process_parameters(parameters):
        """ external parameters and injected beam partameters """
        res = dict()
        """ kinetic energy @ injection """
        W_i   = parameters.get("Win",None) 
        Tk_i  = parameters.get("Tkin",None) 
        res['injection_energy'] = UTIL.PARAMS['injection_energy']  # default is 50 MeV
        if W_i != None:
            res['injection_energy'] = W_i
        elif Tk_i != None:
            res['injection_energy'] = Tk_i

        """ longitudinal energy spread @ injection: {Dphi,w}-space """
        DT2T_i = parameters.get('DT2T',None) 
        DT2T_i = parameters.get('DT2T',None) 
        res['DT2T_i']    = 0.01 # default 1%
        if DT2T_i != None:
            res['DT2T_i'] = DT2T_i
        elif DT2T_i != None:
            res['DT2T_i'] = DT2T_i
        res['Dphi0_i']    = radians(parameters.get('DPHI0',10.)) # default [rad]

        """ transverse beam parameters """
        res['emitx_i']          = parameters.get('emitx_i',1E-6) # [m*rad]
        res['betax_i']          = parameters.get('betax_i',1.)   # [m]
        res['alfax_i']          = parameters.get('alfax_i',0.)
        res['emity_i']          = parameters.get('emity_i',1E-6)
        res['betay_i']          = parameters.get('betay_i',1.)
        res['alfay_i']          = parameters.get('alfay_i',0.)
        # transverse Twiss @ entrance
        res['twiss_x_i']        = UTIL.Twiss(res['betax_i'], res['alfax_i'],res['emitx_i'])
        res['twiss_y_i']        = UTIL.Twiss(res['betay_i'], res['alfay_i'],res['emity_i'])
        UTIL.DEBUG_OFF(f"twiss_x_i {res['twiss_x_i']()}")
        UTIL.DEBUG_OFF(f"twiss_y_i {res['twiss_y_i']()}")
        # initial dispersion @ entrance
        res['dx_i']              = parameters.get('dx_i',0.)
        res['dxp_i']             = parameters.get('dxp_i',0.)
        # supplemental parameters
        res['nbsigma']          = parameters.get('nbsigma',3)
        res['lattice_version']  = parameters.get('lattvers','not given')
        res['thins']            = parameters.get('thins',1)
        res['input_file']       = None
        """ 
        longitudinal beam emittance @ entrance
            Phase space ellipse parameters nach T.Wangler (6.47-48) pp.185
            (w/w0)**2 + (Dphi/Dphi0)**2 = 1
            emitw = w0*Dphi0 = ellipse_area/pi
        """
        Dphi0   = res['Dphi0_i']
        DT2T    = res['DT2T_i']
        T       = res['injection_energy']
        E0      = UTIL.PARAMS['proton_mass']
        w0      = T/E0*DT2T # Wrangler's definition of w (pp.176)
        emit_w  = Dphi0*w0 # emittance in {Dphi,w}-space
        beta_w  = emit_w/w0**2 # twiss-beta in {Dphi,w}-space
        alfa_w  = 0. # immer in {Dphi,w}-space
        res['w0_i']     = w0
        res['emitw_i']  = emit_w
        res['alfaw_i']  = alfa_w
        res['betaw_i']  = beta_w
        # longitudinal TWiss @ entrance in {Dphi,w}-space
        res['twiss_w_i'] = UTIL.Twiss(res['betaw_i'], res['alfaw_i'],res['emitw_i'])
        return res

    def process_elements(elements):
        return elements

    def make_lattice(elementIDs):
        UTIL.DEBUG_OFF(elementIDs)
        lattice = LAT.Lattice(descriptor=UTIL.PARAMS.get('descriptor'))
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
            if isinstance(instance, (ELM.Node)):  # add Node objects only!
                instances.append(instance)
        UTIL.DEBUG_OFF(instances)
        for instance in instances:
            lattice.add_node(instance)
        return lattice   # the complete lattice
    
# >>>>> factory >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    with open(input_file,'r') as fileobject:
        try:
            in_data = yaml.load(fileobject,Loader=yaml.Loader)
        except Exception as ex:
            raise(UserWarning(wrapRED('File inputError: {}'.format(str(ex)))))
            sys.exit(1)
    fileobject.close()
    UTIL.DEBUG_OFF(in_data)

    # call lattice parser, get results
    results = LP2.parse(in_data)
    UTIL.PARAMS['descriptor'] = results.DESCRIPTOR  # get DESCRIPTOR from parsed results

    flags = process_flags(results.FLAGS)
    UTIL.FLAGS.update(flags)
    UTIL.DEBUG_OFF('global FLAGS after process_flags():',UTIL.FLAGS)

    parameters = process_parameters(results.PARAMETERS)
    parameters['input_file'] = input_file
    UTIL.PARAMS.update(parameters)
    UTIL.DEBUG_OFF('global PARAMS after process_parameters():',UTIL.PARAMS)

    elements = process_elements(results.ELEMENTS)
    UTIL.ELEMENTS = elements
    UTIL.DEBUG_OFF('ELEMENTS after process_elements():',UTIL.ELEMENTS)

    lat_elementIDs = results.LAT_ELMIDs
    UTIL.DEBUG_OFF('LAT_ELMIDs after process_elements():',lat_elementIDs)
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
        input_file = "unittests/TT28_base.yml"
        lattice = factory(input_file)
        print(F'\u26dd  FINAL kinetic energy {lattice.seq[-1].ref_track[UTIL.EKOO]:.3f} [MeV] \u26dd')
        UTIL.DEBUG_ON('lattice_generator'   ,lattice.toString())
        UTIL.DEBUG_ON('SUMMARY in factory()',UTIL.SUMMARY)
        UTIL.DEBUG_ON('FLAGS in factory()'  ,UTIL.FLAGS)
        UTIL.DEBUG_ON('PARAMS in factory()' ,UTIL.PARAMS)

if __name__ == '__main__':
    unittest.main()
