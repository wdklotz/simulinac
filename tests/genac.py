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
from math import sqrt,degrees
import logging

from elements import MDIM,XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
import elements as ELM
from setutil import CONF,SUMMARY,Particle,Proton,dictprnt,collect_summary
from lattice_generator import parse_yaml_and_fabric
from tracks import track_soll

def quadF(elm):
    res = '{}'.format('QUADRUPO')
    xl = elm.length*100.                          #XL [cm]
    rg = CONF['quad_bore_radius']                 #RG [m]
    bq =  +10.*elm.k0*elm.particle.brho*rg        #BQ kGauss
    res += '\n{} {} {}'.format(xl,bq,rg*100.)
    return res

def quadD(elm):
    res = '{}'.format('QUADRUPO')
    xl = elm.length*100.                          #XL [cm]
    rg = CONF['quad_bore_radius']                 #RG [m]
    bq = -10.*elm.k0*elm.particle.brho*rg         #BQ kGauss (negative B-Feld!)
    # print('xl rg k0 brho bq',xl,rg,elm.k0,elm.particle.brho,bq)
    res += '\n{} {} {}'.format(xl,bq,rg*100.)
    return res

def drift(elm):
    res = '{}'.format('DRIFT')
    res += '\n{}'.format(elm.length*100.)         #LD
    return res

def cavity(elm):
    res = '{}'.format('CAVSC')
    res += '\n{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
        0.,                               #ETCELL(1)
        0.,                               #ETCELL(2)
        0.,                               #ETCELL(3)
        elm.length*100.,                  #ETCELL(4) cell length [cm]
        0.7,                              #ETCELL(5) T
        0.0007,                           #ETCELL(6) T'
        0.,                               #ETCELL(7) S
        0.,                               #ETCELL(8) S'
        0.,                               #ETCELL(9)
        0.,                               #ETCELL(10)
        elm.u0/CONF['spalt_laenge'],      #ETCELL(11) E0Z [MV/m]
        degrees(elm.phis),                #ETCELL(12) RF phase [deg]
        0.,                               #ETCELL(13)
        0.,                               #ETCELL(14) T''
        elm.freq*1.e-6,                   #ETCELL(15) f [MHz]
        0.01,                             #ETCELL(16) att. factor E0Z
        )
    return res

def rfgap(attributes):
    return 'rfgap not implemented'

#create logger
logger = logging.getLogger("logger")
logger.setLevel(logging.DEBUG)
#create console handler
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
#create formatter
formatter = logging.Formatter("%(levelname)s: %(filename)s[%(lineno)d] %(message)s")
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)

lattice = parse_yaml_and_fabric('fodo_with_10cav_per_RF(3).yml')
Particle.soll = Proton(CONF['injection_energy'])
soll_track = track_soll(lattice)               # track soll Teilchen hier!  (WICHTIG)

dynacIN  =   '{}'.format(CONF['lattice_version'])

dynacIN += '\n{}'.format('GEBEAM')
dynacIN += '\n{} {}'.format(1,1)                             #LAW ITWISS
dynacIN += '\n{}'.format(CONF['frequenz'])                   #FH
dynacIN += '\n{}'.format(9000)                               #IMAX
dynacIN += '\n{} {} {} {} {} {}'.format(0.,0.,0.,0.,0.,0.)   #CENTRE(I)(I=1 to 6)
dynacIN += '\n{} {} {}'.format(0.,0.780,1.)                  #ALPHAX BETAX EMITX
dynacIN += '\n{} {} {}'.format(0.,2.373,1.)                  #ALPHAY BETAY EMITY
dynacIN += '\n{} {} {}'.format(0.0179,0.77283,512.82)        #ALPHAZ BETAZ EMITZ

dynacIN += '\n{}'.format('INPUT')
dynacIN += '\n{} {} {}'.format(CONF['proton_mass'],1.,1.)    #UEM ATM Q
dynacIN += '\n{} {}'.format(CONF['injection_energy'],0.0)    #ENEDEP TOF

# dynacIN += '\n{}'.format('EMIPRT')
# dynacIN += '\n{}'.format(0)                                  #IEMQESG

# dynacIN += '\n{}'.format('SCDYNAC')
# dynacIN += '\n{}'.format(3)                                  #ISCSP
# dynacIN += '\n{} {}'.format(38.,3.)                          #BEAMC SCE10
# dynacIN += '\n{}'.format(0)                                  #RDCF

# limits=dict(xlim1=0.3,ylim1=1.0,xlim2=0.3,ylim2=1.0,xlim3=0.3,ylim3=0.3,xlim4=60.,ylim4=0.06)
limits=dict(xlim1=0.3,ylim1=3.,xlim2=0.3,ylim2=3.,xlim3=0.3,ylim3=0.3,xlim4=40.,ylim4=0.05)
dynacIN += '\n{}'.format('EMITGR')
dynacIN += '\n{}'.format('BEAM AT INPUT')                    #TITLE
dynacIN += '\n{} {}'.format(0,5)                             #IDWDP RMSMTP
dynacIN += '\n{} {} {} {} {} {} {} {}'.format(               #XLIM1 YLIM1 XLIM2 YLIM2 XLIM3 YLIM3 XLIM4 YLIM4
        limits['xlim1'],
        limits['ylim1'],
        limits['xlim2'],
        limits['ylim2'],
        limits['xlim3'],
        limits['ylim3'],
        limits['xlim4'],
        limits['ylim4'])

for cnt,ipos in enumerate(lattice.seq[:2]):
    elm,s0,s1 = ipos                           #element, start, end
    tks = (soll_track.point_at(cnt))[EKOO]     #tk
    ts0 = (soll_track.point_at(cnt))[SKOO]     #s0
#   logger.debug('{}\t(adr {}) label \'{}\'\ttk {:.4f} s0 {:.4f} s1 {:.4f}'.format(elm.__class__,id(elm),elm.label,tks,s0,s1))
    dynacIN += '\n;---{}---'.format(cnt)
    if isinstance(elm,ELM.QD):
        dynacIN += '\n{}'.format(quadD(elm))
    elif isinstance(elm,ELM.QF):
        dynacIN += '\n{}'.format(quadF(elm))
    elif isinstance(elm,ELM.D):
        dynacIN += '\n{}'.format(drift(elm))
    elif isinstance(elm,ELM.RFC):
        dynacIN += '\n{}'.format(cavity(elm))
    elif isinstance(elm,ELM.RFG):
        dynacIN += '\n{}'.format(rfgap(elm))
    elif isinstance(elm,ELM.MRK):
        dynacIN += '\n;MRK(tk[MeV]={}, s[m]={})'.format(tks,s0)
    else:
        pass

limits=dict(xlim1=0.3,ylim1=90.,xlim2=0.3,ylim2=160.,xlim3=0.3,ylim3=0.3,xlim4=60.,ylim4=0.5)
dynacIN += '\n{}'.format('EMITGR')
dynacIN += '\n{}'.format('BEAM AT OUTPUT')                   #TITLE
dynacIN += '\n{} {}'.format(0,5)                             #IDWDP RMSMTP
dynacIN += '\n{} {} {} {} {} {} {} {}'.format(               #XLIM1 YLIM1 XLIM2 YLIM2 XLIM3 YLIM3 XLIM4 YLIM4
        limits['xlim1'],
        limits['ylim1'],
        limits['xlim2'],
        limits['ylim2'],
        limits['xlim3'],
        limits['ylim3'],
        limits['xlim4'],
        limits['ylim4'])

dynacIN += '\n{}'.format('STOP')
print(dynacIN)
