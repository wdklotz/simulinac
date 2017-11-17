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
from math import degrees

from xml_utils.XmlDataAdaptor import XmlDataAdaptor
from setutil import FLAGS,PARAMS,DEBUG
from lattice_generator import factory
from lattice import Lattice
import elements as ELM
import TTFG as TTF


# print(dir())

def generator(dir='yml/', file='25_09_2017_versuche_70_200MeV', ext='.yml', EzFile=None, aperture=None):
    input      = '{}{}{}'.format(dir,file,ext)
    lattice = factory(input)

    root_da  = XmlDataAdaptor(name='Alceli')
    sections = lattice.get_sections()       #sections is a [Section,...]

    for section in sections:
        if len(section.seq) == 0: continue
        sec      = section.get_name()
        sec_da   = root_da.createChild(sec)
        cavs_da  = sec_da.createChild('Cavities')
        gap_cnt  = 0
        quad_cnt = 0
        for node in section.seq:
            if isinstance(node,(ELM.QFth,ELM.QFthx,ELM.QDth,ELM.QDthx)):
                raise RuntimeError('thin QUADs not implememted in pyOrbit. Use thick QUADs only! -- STOP!')
                sys.exit(1)
            elif isinstance(node,(ELM.GAP,ELM.RFC,TTF.TTFG)):
                raise RuntimeError('GAP,RFC,TTFG not compatible with pyOrbit. Use RFG only! -- STOP!')
                sys.exit(1)
            elif isinstance(node,(ELM.SD,ELM.RD,ELM.WD)):
                pass
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
                    accelm_da.setValue('type','QUAD')
                    accelm_da.setValue('name','{}:{}:{}'.format(node.label,sec,quad_cnt))

                    k0 = node.k0
                    if isinstance(node,ELM.QD): k0 = -k0
                    Bgrad = k0*node.particle.brho

                    par_da.setValue('field', Bgrad)
                    par_da.setValue('aperture', aperture)
                    par_da.setValue('aprt_type', 1)

                elif isinstance(node,ELM.RFG):
                    gap_cnt += 1
                    accelm_da.setValue('type','RFGAP')
                    accelm_da.setValue('name','{}:{}:{}'.format(node.label,sec,gap_cnt))
                    ttf_da = accelm_da.createChild('TTFs')

                    phiSoll = degrees(node.phis)
                    E0L = node.u0*1.e-3
                    E0TL = E0L*node.tr

                    par_da.setValue('EOL', E0L)
                    par_da.setValue('EOTL', E0TL)
                    par_da.setValue('EzFile', EzFile)
                    par_da.setValue('aperture', aperture)
                    par_da.setValue('aprt_type', 1)
                    par_da.setValue('cavity', '{}:{}:{}'.format('pillbox',sec,gap_cnt))
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
                    cavity_da.setValue('name', '{}:{}:{}'.format('pillbox',sec,gap_cnt))
                    cavity_da.setValue('pos', sm)
                else:
                    pass

#     DEBUG('root_da.makeXmlText()\n',root_da.makeXmlText())
    root_da.writeToFile('{}{}'.format(file,'.xml'))
    print('----------------------------XmlGenerator for pyOrbit -----')
    print('Result in file ==> {}{}'.format(file,'.xml'))
    print('----------------------------------------------------------')
    return

if __name__ == '__main__':
    EzFile = './SF_WDK2g44.TBL'
    aperture = PARAMS['quad_bore_radius']
    generator(EzFile=EzFile, aperture=aperture)
