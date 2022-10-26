#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
__version__='vv10.22.7'
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
from math import degrees, sqrt
import warnings
# import pickle
import argparse

from xml_utils.XmlDataAdaptor import XmlDataAdaptor
from setutil import FLAGS,PARAMS,DEBUG_ON,DEBUG_OFF
from lattice_generator import factory
from lattice import Lattice
import elements as ELM

def generator(file=None):
    input   = file
    output  = '{}/{}.{}'.format('PyOrbitLattice','generated','xml')
    print("this run ==> input {}, output {}".format(input,output))

    lattice = factory(input)

    # pickle test for whole lattice
    # pickled_lattice   = pickle.dumps(lattice)
    # unpickled_lattice = pickle.loads(pickled_lattice)
    # lattice           = unpickled_lattice
    # DEBUG_ON('lattice length',lattice.length)

    # waccept(lattice.first_gap)
    lattice.first_gap.waccept()

    # select parameters for <PARAMS/>-tag
    parameter = dict(
            clight           = PARAMS['clight'],              # [m/s]
            proton_mass      = PARAMS['proton_mass'],         # MeV
            injection_energy = PARAMS['injection_energy'],    # [MeV]
            lattice_version  = PARAMS['lattice_version'],
            # EzAvg            = PARAMS['EzAvg'],               # [MV/m]
            # gap              = PARAMS['gap'],                 # [m]
            # cavity_laenge    = PARAMS['cavity_laenge'],       # [m]
            # phisoll          = PARAMS['phisoll'],             # [deg]
            # frequenz         = PARAMS['frequenz'],            # [Hz]
            # emitx_i          = PARAMS['emitx_i'],             # [m*rad]
            # emity_i          = PARAMS['emity_i'],             # [m*rad]
            # emitz_i          = PARAMS['emitz'],               # [m*rad]
            # emitw_i          = PARAMS['emitw'],               # [rad]
            # betax_i          = PARAMS['betax_i'],             # [m/rad]
            # betay_i          = PARAMS['betay_i'],             # [m/rad]
            # betaz_i          = PARAMS['betaz'],               # [m/rad]
            # alfax_i          = PARAMS['alfax_i'],             # []
            # alfay_i          = PARAMS['alfay_i'],             # []
            # alfaz_i          = PARAMS['alfaz'],               # []
            # z0               = PARAMS['z0'],                  # [m]
            )
    root_da        = XmlDataAdaptor(name='Alceli')
    params_dict_da = root_da.createChild('PARAMS')
    lattice_da     = root_da.createChild('Lattice')
    cavs_da        = lattice_da.createChild('Cavities')
    for key in parameter:
        value = parameter[key]
        params_dict_da.setValue(key,value)
    lattice_da.setParam('length',lattice.length)

    print('{}'.format(root_da.makeXmlText()))
    print('generating xml for pyOrbit...')

    # iterate through the lattice
    accElement_cnt = 0
    gap_cnt        = 0
    quadF_cnt      = 0
    quadD_cnt      = 0
    quad_cnt       = 0
    drift_cnt      = 0
    nodes  = iter(lattice)
    for node in nodes:
        if isinstance(node,(ELM.GAP, ELM.RFC)):
            raise RuntimeError('GAP,RFC,TTFG not compatible with pyOrbit. Use RFG only! -- STOP!')
            sys.exit(1)
        elif isinstance(node,(ELM.SD, ELM.RD)):
            warnings.warn('dipoles type SD or RD are ignored in generation of pyOrbit lattice')
        elif isinstance(node,(ELM.QF,ELM.QD,ELM.RFG)):
            s0 = node.position[0]  #from
            sm = node.position[1]  #middle position
            s1 = node.position[2]  #to
            accelm_da = lattice_da.createChild('accElement')
            accElement_cnt += 1
            accelm_da.setValue('length',node.length)
            accelm_da.setValue('pos',sm)
            par_da = accelm_da.createChild('parameters')

            if isinstance(node,(ELM.QF,ELM.QD)):
                quad_cnt  += 1
                quadF_cnt += 1
                name = '{}:{}'.format(node.label,quad_cnt)
                # k0 = node.k02
                # Bgrad = sqrt(k0)*node.particle.brho
                Bgrad = node.grad
                if isinstance(node,ELM.QD): 
                    quadF_cnt -= 1
                    quadD_cnt += 1
                    Bgrad = -Bgrad
                aperture = node.aperture

                accelm_da.setValue('type','QUAD')
                accelm_da.setValue('name',name)

                par_da.setValue('field', Bgrad)
                par_da.setValue("dBdr", Bgrad)
                par_da.setValue('aperture', aperture)
                par_da.setValue('aprt_type', 1)
            elif isinstance(node,ELM.RFG):
                gap_cnt += 1
                name = '{}:{}'.format(node.label,gap_cnt)
                ttf_da = accelm_da.createChild('TTFs')

                # phiSoll = degrees(node.phis) + 180. # pyOrbit's soll phase ~135 [deg]!
                phiSoll = degrees(node.phisoll)
                E0L  = node.EzAvg*node.gap*1.e-3   # pyOrbit [Gev]
                E0TL = E0L*node.ttf                # use this when energy is adjusted
                E0TL = E0L*0.8575                  # some reasonable(?) average
                aperture = node.aperture
                fieldtab = node.fieldtab           # file name of SF table
                name = '{}:{}'.format('pillbox',gap_cnt)
                accelm_da.setValue('type','RFGAP')
                accelm_da.setValue('name',name)

                par_da.setValue('E0L', E0L)
                par_da.setValue('E0TL', E0TL)
                par_da.setValue('aperture', aperture)
                par_da.setValue('aprt_type', 1)
                par_da.setValue('cavity', name)
                par_da.setValue('mode', 0)
                par_da.setValue('phase', phiSoll)
                par_da.setValue('fieldtab',fieldtab)

                ttf_da.setValue('beta_max', 1.0)
                ttf_da.setValue('beta_min', 0.0)

                polyT_da  = ttf_da.createChild('polyT')
                polyS_da  = ttf_da.createChild('polyS')
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
                cavity_da.setValue('frequency',node.freq)
                cavity_da.setValue('name', name)
                cavity_da.setValue('pos', sm)
            DEBUG_OFF('accelm_da: {}'.format(accelm_da.makeXmlText()))
        elif isinstance(node,(ELM.D,ELM.DKD)):
            drift_cnt += 1
        else:
            pass
    print("accElement_cnt {}".format(accElement_cnt))
    print("gap_cnt        {}".format(gap_cnt))
    print("quadF_cnt      {}".format(quadF_cnt))
    print("quadD_cnt      {}".format(quadD_cnt))
    print("quad_cnt       {}".format(quad_cnt))
    print("drift_cnt      {}".format(drift_cnt))

    root_da.writeToFile(output)
    return

if __name__ == '__main__':
    # r_limit = sys.getrecursionlimit()
    # DEBUG_OFF(f'recursionlimit {r_limit}')
    # sys.setrecursionlimit(r_limit*60)
    # r_limit = sys.getrecursionlimit()
    # DEBUG_ON(f'recursionlimit {r_limit}')

    # use ArgumentParser to put result in 'args'
    parser = argparse.ArgumentParser()
    group  = parser.add_mutually_exclusive_group()
    group.add_argument ("--file", default="yml/simuINwork.yml",   help="lattice input-file:default simuINwork.yml")
    args = vars(parser.parse_args())
    DEBUG_OFF(args)

    file_name = args['file']
    generator(file = file_name)
