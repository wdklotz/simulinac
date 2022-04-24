#!/Users/klotz/anaconda3/bin/python3.6
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
import os
import sys
from math import degrees
import warnings

from xml_utils.XmlDataAdaptor import XmlDataAdaptor
from setutil import FLAGS,PARAMS,DEBUG,waccept
from lattice_generator import factory
from lattice import Lattice
import elements as ELM

# DEBUG
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass

DEBUG_GEN    = DEBUG_OFF
# DEBUG_GEN(lineno(),dir())

def generator(dir='yml', file='simuIN', ext='yml', EzFile=None):
    input   = '{}/{}.{}'.format(dir,file,ext)
    output  = '{}/{}.{}'.format('.','lattice','xml')

    # lattice and longitudinal paramters at entrance
    lattice = factory(input)
    waccept(lattice.first_gap)

    root_da   = XmlDataAdaptor(name='Alceli')
    sections  = lattice.get_sections() #sections is [section-1,...,section-N]

    # transfer a selection of parameters to <PARAMS/>-tag
    parameter = dict(
            clight           = PARAMS['clight'],# [m/s]
            proton_mass      = PARAMS['proton_mass'],         # MeV
            # EzAvg            = PARAMS['EzAvg'],               # [MV/m]
            gap              = PARAMS['gap'],                 # [m]
            cavity_laenge    = PARAMS['cavity_laenge'],       # [m]
            # phisoll          = PARAMS['phisoll'],             # [deg]
            frequenz         = PARAMS['frequenz'],            # [Hz]
            injection_energy = PARAMS['injection_energy'],    # [MeV]
            emitx_i          = PARAMS['emitx_i'],             # [m*rad]
            emity_i          = PARAMS['emity_i'],             # [m*rad]
            emitz_i          = PARAMS['emitz'],               # [m*rad]
            emitw_i           = PARAMS['emitw'],              # [rad]
            betax_i          = PARAMS['betax_i'],             # [m/rad]
            betay_i          = PARAMS['betay_i'],             # [m/rad]
            betaz_i          = PARAMS['betaz'],               # [m/rad]
            alfax_i          = PARAMS['alfax_i'],             # []
            alfay_i          = PARAMS['alfay_i'],             # []
            alfaz_i          = PARAMS['alfaz'],               # []
            lattice_version  = PARAMS['lattice_version'],
            z0               = PARAMS['z0'],                  # [m]
                     )
    params_dict_da = root_da.createChild('PARAMS')
    for key in parameter:
        value = parameter[key]
        params_dict_da.setValue(key,value)

    print('-----XmlGenerator for pyOrbit -----\n{}'.format(root_da.makeXmlText()))
    print('wait...')

    # loop over sections
    for section in sections:
        if len(section.seq) == 0: continue
        sec      = "S"+section.name
        sec_da   = root_da.createChild(sec)
        sec_da.setValue('length',section.length)
        sec_da.setValue('name',sec)
        DEBUG_GEN('sec_da: {}'.format(sec_da.makeXmlText()))
        cavs_da  = sec_da.createChild('Cavities')
        gap_cnt  = 0
        quad_cnt = 0

        for node in section.seq:
            if isinstance(node,(ELM.QFth, ELM.QFthx, ELM.QDth ,ELM.QDthx)):
                raise RuntimeError('thin QUADs not implememted in pyOrbit. Use thick QUADs only! -- STOP!')
                sys.exit(1)
            elif isinstance(node,(ELM.GAP, ELM.RFC)):
                raise RuntimeError('GAP,RFC,TTFG not compatible with pyOrbit. Use RFG only! -- STOP!')
                sys.exit(1)
            elif isinstance(node,(ELM.SD, ELM.RD)):
                warnings.warn('SD,RD are ignored in generation of pyOrbit lattice')
            elif isinstance(node,(ELM.QF,ELM.QD,ELM.RFG)):
                s0 = node.position[0]  #from
                sm = node.position[1]  #middle position
                s1 = node.position[2]  #to
                accelm_da = sec_da.createChild('accElement')
                accelm_da.setValue('length',node.length)
                accelm_da.setValue('pos',sm)
                par_da = accelm_da.createChild('parameters')

                if isinstance(node,(ELM.QF,ELM.QD)):
                    quad_cnt += 1
                    name = '{}:{}'.format(node.label,quad_cnt)
                    k0 = node.k0
                    if isinstance(node,ELM.QD): k0 = -k0
                    Bgrad = k0*node.particle.brho
                    aperture = node.aperture

                    accelm_da.setValue('type','QUAD')
                    accelm_da.setValue('name',name)

                    par_da.setValue('field', Bgrad)
                    par_da.setValue('aperture', aperture)
                    par_da.setValue('aprt_type', 1)

                elif isinstance(node,ELM.RFG):
                    gap_cnt += 1
                    name = '{}:{}'.format(node.label,gap_cnt)
                    ttf_da = accelm_da.createChild('TTFs')

                    # phiSoll = degrees(node.phis) + 180. # pyOrbit's soll phase ~135 [deg]!
                    phiSoll = degrees(node.phis)
                    E0L  = node.EzAvg*node.gap*1.e-3   # pyOrbit [Gev]
                    E0TL = E0L*node.ttf                # use this when energy is adjusted
                    E0TL = E0L*0.8575                  # some reasonable(?) average
                    aperture = node.aperture
                    name = '{}:{}'.format('pillbox',gap_cnt)
                    accelm_da.setValue('type','RFGAP')
                    accelm_da.setValue('name',name)

                    par_da.setValue('E0L', E0L)
                    par_da.setValue('E0TL', E0TL)
                    par_da.setValue('EzFile', EzFile)
                    par_da.setValue('aperture', aperture)
                    par_da.setValue('aprt_type', 1)
                    par_da.setValue('cavity', name)
                    par_da.setValue('mode', 0)
                    par_da.setValue('phase', phiSoll)

                    ttf_da.setValue('beta_max', 1.0)
                    ttf_da.setValue('beta_min', 0.0)

                    polyT_da = ttf_da.createChild('polyT')
                    polyS_da = ttf_da.createChild('polyS')
                    polyTP_da = ttf_da.createChild('polyTP')
                    polySP_da = ttf_da.createChild('polySP')
                    polyT_da.setValue('order', 2)
                    polyS_da.setValue('order', 2)
                    polyTP_da.setValue('order', 2)
                    polySP_da.setValue('order', 2)
                    polyT_da.setValue('pcoefs', (1.,0.,0.))
                    polyS_da.setValue('pcoefs', (1.,0.,0.))
                    polyTP_da.setValue('pcoefs', (1.,0.,0.))
                    polySP_da.setValue('pcoefs', (1.,0.,0.))

                    cavity_da = cavs_da.createChild('Cavity')
                    cavity_da.setValue('ampl', 1.)
                    cavity_da.setValue('frequency', PARAMS['frequenz'])
                    cavity_da.setValue('name', name)
                    cavity_da.setValue('pos', sm)
                else:
                    pass
            DEBUG_GEN('accelm_da: {}'.format(accelm_da.makeXmlText()))

    root_da.writeToFile(output)
    print('Input from file ==> {}'.format(input))
    print('Result in  file ==> {}'.format(output))
    return

if __name__ == '__main__':
    EzFile = './SF_WDK2g44.TBL'
    generator(EzFile = EzFile)
