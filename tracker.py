#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
___version___='v7.0.4a2'
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
import numpy as NP
import matplotlib.pyplot as plt
import time
from string import Template
from math import sqrt, degrees, radians

from lattice_generator import factory
import elements as ELM
import marker_actions as MRK
from setutil import DEBUG, PARAMS, FLAGS, dictprnt, sigmas, K, PARAMS, waccept
from setutil import WConverter
from bunch import BunchFactory, Gauss1D, Track, Tpoint, Bunch
from trackPlot import poincarePlot

# DEBUGGING
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass

def scatterPlot(live_lost, abszisse, ordinate, text, minmax=(1.,1.)):
    """ 
    Prepare the Poincaré section plot 
    """
    live_bunch, lost_bunch = live_lost
    txt = ('{} initial'.format(text),'{} final'.format(text))
    initial = 0; final   = -1
    x=[]; y=[]; xlost=[]; ylost=[]
    loc = initial                   # INITIAL
    fig = plt.figure(loc)
    nbprt = live_bunch.nbofparticles()+lost_bunch.nbofparticles()
    box = '{} {} particles'.format(txt[loc],nbprt)
    for particle in iter(live_bunch): # particles
        track  = particle.track
        tpoint = track.getpoints()[loc]
        point  = tpoint()
        x.append(point[abszisse])
        y.append(point[ordinate])
    for particle in iter(lost_bunch): # lost particles
        track  = particle.track
        tpoint = track.getpoints()[loc]
        point  = tpoint()
        xlost.append(point[abszisse])
        ylost.append(point[ordinate])
    xmax = max(x)
    xmin = min(x)
    ymax = max(y)
    ymin = min(y)
    minmax = (xmax-xmin,ymax-ymin)
    # plt.scatter(x,y,s=1)   # bare scatter plot
    poincarePlot((x,y),(xlost, ylost), box, max = minmax, projections = (1,1))

    x=[]; y=[]
    loc = final                     # FINAL
    fig = plt.figure(loc)
    nbprt = live_bunch.nbofparticles()
    box = '{} {} particles'.format(txt[loc],nbprt)
    for particle in iter(live_bunch): # particles
        track  = particle.track
        tpoint = track.getpoints()[loc]
        point  = tpoint()
        x.append(point[abszisse])
        y.append(point[ordinate])
    xmax = max(x)
    xmin = min(x)
    ymax = max(y)
    ymin = min(y)
    minmax = (xmax-xmin,ymax-ymin)
    # plt.scatter(x,y,s=1)   # bare scatter plot
    poincarePlot((x,y),(0,0), box, max = minmax, projections = (1,1))
    return
    
def projection(live_lost, show, abszisse= K.z, ordinate= K.zp):
    """ 
    2D phase space projections of Poincaré sections 
    """
    symbols = ("x","x'","y","y'","z","$\Delta$p/p")
    text    = '{}-{}'.format(symbols[abszisse],symbols[ordinate])
    fig     = scatterPlot(live_lost, abszisse=abszisse, ordinate=ordinate, text=text)
    if show: plt.show()

def frames(lattice, save, skip):
    if save:
        nscnt = 0
        for node in iter(lattice):
            if isinstance(node,MRK.PoincareAction):
                nscnt += 1
                if nscnt%skip == 0:
                    node.do_action(nscnt)
    
def progress(tx):
    template = Template('$tx1 $tx2 $tx3 $tx4')
    res = template.substitute(tx1=tx[0] , tx2=tx[1] , tx3=tx[2] , tx4=tx[3] )
    print('{}\r'.format(res),end='')

#todo: aperture check
#todo: useaper flag
def track_node(node,particle):
    """
    Tracks a particle through a node
    """
    track   = particle.track
    last_tp = track.getpoints()[-1]
    try:
        new_point = node.map(last_tp())
        new_tp    = Tpoint(point=new_point)
    except ValueError as ex:
        lost = True
        track.removepoint(last_tp)
        return lost
    # check Dp2p-acceptance
    if abs(new_point[K.zp]) < PARAMS['Dp2pAcceptance']:
        lost = False
    else:
        lost = True
    # check z-acceptance
    if abs(new_point[K.z]) < PARAMS['zAcceptance']:
        lost = False
    else:
        lost = True
    if not lost:
        if track.nbofpoints() > 1:
            # !!DISCARD!! last point
            track.removepoint(last_tp)
        track.addpoint(new_tp)
        if isinstance(node,MRK.PoincareAction) and node.has_action('scatter'):
            node.add_track_point(new_tp)
    particle.lost = lost
    return lost
    
def track(lattice,bunch):
    """
    Tracks a bunch of particles through the lattice using maps
    - lattice is a list of elements (class _Node)
    - bunch (class Bunch) is a list of particles (class Particle)
    - each particle in a bunch has a track=paticle['track']
    - track (class Track) is a list of points (class Tpoint)
    
    Input: lattice , bunch
    """
    zeuge          = ('*\u007C*','**\u007C','*\u007C*','\u007C**')  # *|*
    tx4            = ' {}% done {}/{} lost/initial'.format(0,0,bunch.nbofparticles())

    nscnt  = 0
    ndcnt  = 0
    lnode  = len(lattice.seq)
    lmod   = int(lnode*0.05)
    nlost  = 0
    nbpart = bunch.nbofparticles()
    lbunch = Bunch()    # lost particles go into this bunch
    for node in iter(lattice):              # nodes
        ndcnt +=1
        for particle in iter(bunch):        # particles
            lost = track_node(node,particle)
            if lost:
                lbunch.addparticle(particle)
                bunch.removeparticle(particle)
                nlost += 1
            # showing track-loop progress
            if (ndcnt+1)%lmod == 0:
                tx4    = ' {}% done {}/{} lost/initial'.format(int((ndcnt/lnode*100.)+1.), nlost, nbpart)
                # tx = ('(track design)', '(track bunch)', zeuge[ndcnt%4], tx4)
                tx = ('(track design)', '(track bunch)', '', tx4)
                progress(tx)
                
    live = nbpart - lbunch.nbofparticles()
    print('TRACKING DONE (live particles {}, lost particles {})               '.format(live,nlost))
    return (bunch,lbunch)

def track_soll(lattice):
    """
    Tracks the reference particle through the lattice and redefines the lattice 
    element parameters according to the energy of the accelerated reference particle.
    """
    soll_track = Track()
    tpoint = Tpoint(NP.array([ 0., 0., 0., 0., 0., 0., PARAMS['sollteilchen'].tkin, 1., 0., 1.]))
    soll_track.addpoint(tpoint)
    for node in iter(lattice):
        pi = soll_track.getpoints()[-1]   # Tpoint at entrance
        """ energy adjustment """
        node.adjust_energy(pi()[K.T])
        """ mapping with soll map """
        pf = node.soll_map(pi())
        tpoint = Tpoint(pf)                  # n+1,Tpoint at exit
        soll_track.addpoint(tpoint)
    return soll_track

def tracker(options):
    """ prepare and launch tracking """
    npart    = options['particles_per_bunch']
    print('-----------------------track_bunch with {} particles---'.format(npart))

    # !!FIRST!! make lattice
    t0       = time.clock()
    filepath = options['input_file']
    lattice  = factory(filepath)

    # No tracking without acceleration
    if FLAGS['dWf'] == 0:
        print('{}'.format('IMPOSSIBLE: no tracking without acceleration?'))
        sys.exit()

    # calculate twiss paramters at entrance
    waccept(lattice.first_gap)
    tkin     = PARAMS['sollteilchen'].tkin
    conv     = WConverter(tkin)
    t1       = time.clock()

    # pull more options
    show     = options['show']
    save     = options['save']
    skip     = options['skip']

    # bunch-configuration from PARAMS
    # {x(x)xp}   standard units
    twx = PARAMS['twiss_x_i']
    betax_i,alfax_i,gammax_i,emitx_i = twx()
    sigma_x   = twx.sigmaH()
    sigma_xp  = twx.sigmaV()
    # {y(x)yp}
    twy = PARAMS['twiss_y_i']
    betay_i,alfay_i,gammay_i,emity_i = twy()
    sigma_y   = twy.sigmaH()
    sigma_yp  = twy.sigmaV()
    # {z(x)Dp2p}  T3D units
    twz = PARAMS['twiss_z_i']
    betaz,alfaz,gammaz,emitz = twz()
    sigma_z    = twz.sigmaH()
    sigma_Dp2p = twz.sigmaV()
    Dp2pmx     = PARAMS['Dp2pmx']
    Dp2p0      = PARAMS['Dp2p0']
    # {Dphi(x)w}  T.Wangler units
    tww = PARAMS['twiss_w_i']
    betaw,alfaw,gammaw,emitw = tww()
    sigma_Dphi  = tww.sigmaH()
    sigma_w     = tww.sigmaV()
    wmx         = PARAMS['wmx']

    # gather for print
    tracker_log = {}
    tracker_log['tkin.......[MeV]'] = tkin
    tracker_log["sigma(x,x')_i.....([m,rad])"] = (sigma_x,sigma_xp)
    tracker_log["sigma(y,y')_i.....([m,rad])"] = (sigma_y,sigma_yp)
    tracker_log["sigma(Dphi,w)_i....([rad,])"] = (sigma_Dphi,sigma_w)
    tracker_log["sigma(z,Dp2p)_i......([m,])"] = (sigma_z,sigma_Dp2p)
    tracker_log['betax_i......[m]'] = betax_i
    tracker_log['betay_i......[m]'] = betay_i
    tracker_log['betaw_i....[rad]'] = betaw
    tracker_log['betaz_i..[m/rad]'] = betaz
    tracker_log['emitx_i......[m]'] = emitx_i
    tracker_log['emity_i......[m]'] = emity_i
    tracker_log['emitw_i....[rad]'] = (emitw, wmx)
    tracker_log['emitz_i......[m]'] = emitz
    tracker_log['Dp2p.........[%]'] = (Dp2p0*1.e2,Dp2pmx*1.e2)
    tracker_log['Dp2p-accptance....[%]'] = PARAMS['Dp2pAcceptance']*1.e2
    tracker_log['z-accpetance.....[mm]'] = PARAMS['zAcceptance']*1.e3
    tracker_log['lattice version......'] = PARAMS['lattice_version']
    dictprnt(tracker_log,'Tracker Log'); print()

    # bunch factory
    bunchfactory = BunchFactory()
    bunchfactory.setDistribution(Gauss1D)
    bunchfactory.setTwiss((twx,twy,twz))
    bunchfactory.setMask(NP.array((1,1,1,1,1,1)))
    bunchfactory.setNumberOfParticles(npart)
    bunchfactory.setReferenceEnergy(PARAMS['injection_energy'])
    bunch = bunchfactory()

    # launch tracking and show final with time
    progress(('(track design)', '', '', ''))
    t2 = time.clock()
    track_soll(lattice)  # <----- track soll
    t3 = time.clock()
    progress(('(track design)', '(track bunch)', '', ''))
    live_lost = track(lattice,bunch) # <----- track bunch returns (live,lost)-bunch
    t4 = time.clock()

    # make 2D projections
    print('FILL PLOTS')
    projection(live_lost, show,  abszisse = K.z, ordinate = K.zp)
    t5 = time.clock()
    print('SAVE FRAMES')
    frames(lattice, save, skip)
    t6 = time.clock()
    # finish up
    print()
    print('total time     >> {:6.3f} [sec]'.format(t5-t0))
    print('parse lattice  >> {:6.3f} [sec] {:4.1f} [%]'.format((t1-t0),(t1-t0)/(t5-t0)*1.e2))
    print('generate bunch >> {:6.3f} [sec] {:4.1f} [%]'.format((t2-t1),(t2-t1)/(t5-t0)*1.e2))
    print('track design   >> {:6.3f} [sec] {:4.1f} [%]'.format((t3-t2),(t3-t2)/(t5-t0)*1.e2))
    print('track bunch    >> {:6.3f} [sec] {:4.1f} [%]'.format((t4-t3),(t4-t3)/(t5-t0)*1.e2))
    print('fill plots     >> {:6.3f} [sec] {:4.1f} [%]'.format((t5-t4),(t5-t4)/(t5-t0)*1.e2))
    print('save frames    >> {:6.3f} [sec] {:4.1f} [%]'.format((t6-t5),(t6-t5)/(t6-t0)*1.e2))

def test0(filepath):
    print('-----------------------------------------Test0---')
    lattice = factory(filepath)
    sollTrack = track_soll(lattice)
    table = sollTrack.as_table()
    DEBUG_TEST0('sollTrack:\n'+table)
    first = sollTrack[0]
    last  = sollTrack[-1]
    DEBUG_TEST0('sollTrack:\n(first): {}\n (last): {}'.format(first.as_str(),last.as_str()))

if __name__ == '__main__':
    DEBUG_TRACK       = DEBUG_OFF
    DEBUG_SOLL_TRACK  = DEBUG_OFF
    DEBUG_TEST0       = DEBUG_ON

    # test0('yml/trackIN.yml')

    print('tracker.py {} on python {}.{}.{}'.format(___version___,sys.version_info.major,sys.version_info.minor,sys.version_info.micro))
    
    # preset files for launch with  m4
    template_file = 'yml/tmpl.yml'          # def.template file
    macros_file   = 'yml/macros.sh'         # def.macro definitions
    input_file    = 'yml/trackIN.yml'       # def.input file

    if sys.platform   == 'win32':
        if len(sys.argv) == 2:
            input_file    = sys.argv[1]
    elif sys.platform == 'darwin' or sys.platform.startswith('linux'):
        if len(sys.argv) == 2:
            input_file    = sys.argv[1]
        else:
            # launch m4
            command = "chmod +x yml/macros.sh"
            command = "{};{} {} > {}".format(command,macros_file,template_file, input_file)
            print('m4 script="{}" template="{}" input="{}"'.format(macros_file,template_file,input_file))
            os.system(command)
    else:
        print('wrong platform')
        sys.exit(1)

    options = {}
    options['input_file']          = input_file
    options['particles_per_bunch'] = 3000
    options['show']                = True
    options['save']                = True
    options['skip']                = 1

    # start the run
    tracker(options)
