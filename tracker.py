#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
___version___='v8.0.0a1'
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
#todo: statistical analalysis of bunch: position, size, etc ...
#todo: how to get the hokey stick?
#todo: check w-acceptance at each node entrance
#todo: no phase damping - why?
import sys,os
import numpy as NP
import matplotlib.pyplot as plt
import time
from string import Template
from math import sqrt, degrees, radians

from lattice_generator import factory
import elements as ELM
import marker_actions as MRK
from setutil import DEBUG,DEBUG_ON,DEBUG_OFF, PARAMS, FLAGS, dictprnt, sigmas, Ktp, PARAMS, waccept
from setutil import WConverter, Functions, TmStamp
from bunch import BunchFactory, Gauss1D, Track, Tpoint, Bunch
# from trackPlot import poincarePlot

def scatterPlot(live_lost, abszisse, ordinate, text, minmax=(1.,1.)):
    """ 
    Plot the scatter plots 
    """
    live_bunch, lost_bunch = live_lost
    txt = ('IN {}'.format(text),'OUT {}'.format(text))
    initial = 0; final   = -1
    # INITIAL DATA
    x=[]; y=[]; xlost=[]; ylost=[]
    loc = initial
    nbprt = live_bunch.nbparticles()+lost_bunch.nbparticles()
    for particle in iter(live_bunch): # particles
        track  = particle.track
        tpoint = track.getpoints()[loc]
        point  = tpoint()
        x.append(point[abszisse]*1e3)
        y.append(point[ordinate]*1e3)
    for particle in iter(lost_bunch): # lost particles
        track  = particle.track
        tpoint = track.getpoints()[loc]
        point  = tpoint()
        xlost.append(point[abszisse]*1e3)
        ylost.append(point[ordinate]*1e3)
    xmax = max([abs(i) for i in x])*1.5
    ymax = max([abs(i) for i in y])*1.5
    # figure with mapping box
    width = 12; height = 6.
    fig   = plt.figure(figsize=(width,height))
    box   = '{} {} particles'.format(txt[loc],nbprt)
    ax    = plt.subplot(121)
    ax.set_title(box)
    # mapping box
    ax.text(0.01, 1.1, PARAMS['mapping'], transform= ax.transAxes, fontsize= 8, bbox= dict(boxstyle='round',facecolor='wheat',alpha=0.5), verticalalignment= 'top')
    plt.xlabel("$10^{-3}$")
    plt.ylabel("$10^{-3}$")
    plt.xlim([-xmax,xmax])
    plt.ylim([-ymax,ymax])
    plt.scatter(x,y,s=1)
    plt.scatter(xlost,ylost,s=1,color='red')
    # poincarePlot((x,y),(xlost, ylost), box, max = minmax, projections = (1,1))

    # FINAL DATA
    x=[]; y=[]
    loc   = final
    nbprt = live_bunch.nbparticles()
    for particle in iter(live_bunch): # particles
        track  = particle.track
        tpoint = track.getpoints()[loc]
        point  = tpoint()
        x.append(point[abszisse]*1e3)
        y.append(point[ordinate]*1e3)
    xmax = max([abs(i) for i in x])*1.5
    ymax = max([abs(i) for i in y])*1.5
    # figure
    box = '{} {} particles'.format(txt[loc],nbprt)
    ax = plt.subplot(122)
    ax.set_title(box)
    plt.xlabel("$10^{-3}$")
    plt.ylabel("$10^{-3}$")
    plt.xlim([-xmax,xmax])
    plt.ylim([-ymax,ymax])
    plt.scatter(x,y,s=1)
    # poincarePlot((x,y),(0,0), box, max = minmax, projections = (1,1))
    return
    
def projections(live_lost):
    """ 
    2D phase space projections IN and OUT
    """
    symbols = ("x","x'","y","y'","z","$\Delta$p/p")
    # (x,xp)
    abszisse = Ktp.x
    ordinate = Ktp.xp
    text    = '{}-{}'.format(symbols[abszisse],symbols[ordinate])
    scatterPlot(live_lost, abszisse=abszisse, ordinate=ordinate, text=text)
    # (y,yp)
    abszisse = Ktp.y
    ordinate = Ktp.yp
    text    = '{}-{}'.format(symbols[abszisse],symbols[ordinate])
    scatterPlot(live_lost, abszisse=abszisse, ordinate=ordinate, text=text)
    # (x,y)
    abszisse = Ktp.x
    ordinate = Ktp.y
    text    = '{}-{}'.format(symbols[abszisse],symbols[ordinate])
    scatterPlot(live_lost, abszisse=abszisse, ordinate=ordinate, text=text)
    # (z,zp)
    abszisse = Ktp.z
    ordinate = Ktp.zp
    text    = '{}-{}'.format(symbols[abszisse],symbols[ordinate])
    scatterPlot(live_lost, abszisse=abszisse, ordinate=ordinate, text=text)
    plt.show()

def frames(lattice, skip):
    """
    2D phase space projection at marker position
    """
    plt.figure()    # new figure instance
    nscnt = 0
    scatter_mrkr = []
    # gather and count markers
    for node in iter(lattice):
        if isinstance(node,MRK.PoincareAction):
            nscnt += 1
            if nscnt%skip == 0:
                scatter_mrkr.append((nscnt,node))
    # make an estimate for x- and y-axis
    d,first_mrkr = scatter_mrkr[0]
    x = [abs(tp()[first_mrkr.xaxis]) for tp in first_mrkr.tpoints]
    y = [abs(tp()[first_mrkr.yaxis]) for tp in first_mrkr.tpoints]
    xmax = max(x)*1.5
    ymax = max(y)*1.5
    # invoke the marker's action
    for nscnt,node in iter(scatter_mrkr):
        node.do_action(nscnt,xmax,ymax)

def loss_plot(lattice,live_lost):
    """
    Plot losses along the lattice
    """
    # figure
    width = 15
    fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)
    fig.suptitle('losses')
    fig.set_size_inches(width,width/2.)
    ax1.set_xlabel('s [m]')
    ax1.set_ylabel('x [mm]')
    ax2.set_xlabel('s [m]')
    ax2.set_ylabel('y [mm]')
    # data
    live_bunch, lost_bunch = live_lost
    for particle in iter(lost_bunch):
        track  = particle.track
        lost = Functions(('s','x','y'))
        for tpoint in iter(track):
            point  = tpoint()
            xlost = point[Ktp.x]*1.e3
            ylost = point[Ktp.y]*1.e3
            s = point[Ktp.S]
            lost.append(s,(xlost,ylost))
        s    = [lost(i,'s') for i in range(lost.nbpoints)]
        ordx = [lost(i,'x') for i in range(lost.nbpoints)]
        ordy = [lost(i,'y') for i in range(lost.nbpoints)]
        ax1.plot(s,ordx)
        ax2.plot(s,ordy)
    lat_plot,d   = lattice.lattice_plot_functions()
    vis_abszisse = [lat_plot(i,'s')     for i in range(lat_plot.nbpoints)]
    vis_ordinate = [lat_plot(i,'viseo') for i in range(lat_plot.nbpoints)]
    vis_zero     = [0.                  for i in range(lat_plot.nbpoints)] # zero line
    ax1.plot(vis_abszisse,vis_ordinate,color='gray')
    ax1.plot(vis_abszisse,vis_zero,color='gray')
    ax2.plot(vis_abszisse,vis_ordinate,color='gray')
    ax2.plot(vis_abszisse,vis_zero,color='gray')
    plt.show()
    return
    
def progress(tx):
    """
    Show progress
    """
    template = Template('$tx1 $tx2 $tx3 $tx4')
    res = template.substitute(tx1=tx[0] , tx2=tx[1] , tx3=tx[2] , tx4=tx[3] )
    # print('\r{}\r'.format(res),end="")
    sys.stdout.write('{}\r'.format(res))

# max track amplitudes set here!
xlim_max  = ylim_max  =  10.e-3
xplim_max = yplim_max =  10.e-3
zlim_max  = zplim_max = 100.e-3
limit     = \
    sqrt(xlim_max**2+xplim_max**2+ylim_max**2+yplim_max**2+zlim_max**2+zplim_max**2)

def track_node(node,particle,options):
    """
    Tracks a particle through a node
    """
    def norm(track):
        track_norm = \
        sqrt(track[0]**2+track[1]**2+track[2]**2+track[3]**2+track[4]**2+track[5]**2)
        return track_norm > limit
    
    track   = particle.track
    last_tp = track.getpoints()[-1]
    lost    = False
    try:
        new_point = node.map(last_tp())
        new_tp    = Tpoint(point=new_point)
    except (ValueError,OverflowError) as ex:
        lost = True
        # track.removepoint(last_tp)   done in track()?
        particle.lost = lost
        return lost

    # limit to stay physicaly reasonable 
    if not FLAGS['useaper']:
        if norm(last_tp()):
            lost = True
            # track.removepoint(last_tp)   done in track()?
            particle.lost = lost
            return

    # check Dp2p-acceptance
    if FLAGS['useaper']:
        if abs(new_point[Ktp.zp]) < PARAMS['Dp2pAcceptance']:
            lost = False
        else:
            lost = True
        # check z-acceptance
        if abs(new_point[Ktp.z]) < PARAMS['zAcceptance']:
            lost = False
        else:
            lost = True

    # check aperture
    if FLAGS['useaper'] and node.aperture != None:
        if abs(new_point[Ktp.x]) < node.aperture and abs(new_point[Ktp.y]) < node.aperture:
            lost = False
        else:
            lost = True

    # if we look for losses we keep all track points
    if not lost:
        if track.nbpoints() > 1 and not options['losses']:
            # !!DISCARD!! last point
            track.removepoint(last_tp)
        track.addpoint(new_tp)
        if isinstance(node,MRK.PoincareAction) and node.has_action('scatter'):
            node.add_track_point(new_tp)
    particle.lost = lost
    return lost
    
def track(lattice,bunch,options):
    """
    Tracks a bunch of particles through the lattice using maps
    - lattice is a list of elements (class _Node)
    - bunch (class Bunch) is a list of particles (class Particle)
    - each particle in a bunch has a track=particle['track']
    - track (class Track) is a list of points (class Tpoint)
    
    Input: lattice , bunch
    """
    zeuge          = ('*\u007C*','**\u007C','*\u007C*','\u007C**')  # *|*
    tx4            = ' {}% done {}/{} lost/initial'.format(0,0,bunch.nbparticles())

    ndcnt  = 0
    lnode  = len(lattice.seq)
    lmod   = int(lnode*0.05)
    nlost  = 0
    nbpart = bunch.nbparticles()
    lbunch = Bunch()    # lost particles go into this bunch
    for node in iter(lattice):              # nodes
        ndcnt +=1
        for particle in iter(bunch):        # particles
            lost = track_node(node,particle,options)
            if lost:
                lbunch.addparticle(particle)
                bunch.removeparticle(particle)
                nlost += 1
            # showing track-loop progress
            if (ndcnt+1)%lmod == 0:
                tx4 = ' {}% done {}% lost'.format(int((ndcnt/lnode*100.)+1.), int(nlost/nbpart*100.))
                # tx = ('(track design)', '(track bunch)', zeuge[ndcnt%4], tx4)
                tx = ('(track design)', '(track bunch)', '', tx4)
                progress(tx)
                
    live = nbpart - lbunch.nbparticles()
    print('\nTRACKING DONE (live particles {}, lost particles {})               '.format(live,nlost))
    return (bunch,lbunch)

def track_soll(lattice):
    """
    Track the reference particle through the lattice and redefines the lattice 
    element parameters according to the energy of the accelerated reference particle.
    """
    soll_track = Track()
    tpoint = Tpoint(NP.array([ 0., 0., 0., 0., 0., 0., PARAMS['sollteilchen'].tkin, 1., 0., 1.]))
    soll_track.addpoint(tpoint)
    for node in iter(lattice):
        pi = soll_track.getpoints()[-1]   # Tpoint at entrance
        """ energy adjustment """
        node.adjust_energy(pi()[Ktp.T])
        """ mapping with soll map """
        pf = node.soll_map(pi())
        tpoint = Tpoint(pf)               # Tpoint at exit
        soll_track.addpoint(tpoint)
    return soll_track

def tracker(options):
    """ 
    Prepare and launch tracking 
    """
    npart    = options['particles_per_bunch']
    print('-----------------------track_bunch with {} particles---'.format(npart))

    # !!FIRST!! make lattice
    t0       = time.process_time()
    filepath = options['input_file']
    lattice  = factory(filepath)

    # No tracking without acceleration
    if FLAGS['dWf'] == 0:
        print('{}'.format('IMPOSSIBLE: no tracking without acceleration!'))
        sys.exit()

    # calculate twiss paramters at entrance
    waccept(lattice.first_gap)
    tkin     = PARAMS['sollteilchen'].tkin
    conv     = WConverter(tkin,lattice.first_gap.freq)
    t1       = time.process_time()

    # pull more options
    show     = options['show']
    save     = options['save']
    skip     = options['skip']
    losses   = options['losses']
    
    # manipulate options
    if losses :
        show = False
        save = False
    if show or save:
        losses = False
    options['show']   = show
    options['save']   = save
    options['losses'] = losses

    # bunch-configuration from PARAMS
    # {x,xp}   standard units
    twx = PARAMS['twiss_x_i']
    betax_i,alfax_i,gammax_i,emitx_i = twx()
    sigma_x   = twx.sigmaH()
    sigma_xp  = twx.sigmaV()
    # {y,yp}
    twy = PARAMS['twiss_y_i']
    betay_i,alfay_i,gammay_i,emity_i = twy()
    sigma_y   = twy.sigmaH()
    sigma_yp  = twy.sigmaV()
    # {z,Dp2p}  T3D units
    twz = PARAMS['twiss_z_i']
    betaz,alfaz,gammaz,emitz = twz()
    sigma_z    = twz.sigmaH()
    sigma_Dp2p = twz.sigmaV()
    Dp2pmx     = PARAMS['Dp2pmx']
    Dp2p0      = PARAMS['Dp2p0']
    # {Dphi,w}  T.Wangler units
    tww = PARAMS['twiss_w_i']
    betaw,alfaw,gammaw,emitw = tww()
    sigma_Dphi  = tww.sigmaH()
    sigma_w     = tww.sigmaV()
    wmx         = PARAMS['wmx']

    # gather for print
    tracker_log = {}
    tracker_log['T-kin.......[MeV]']          = tkin
    tracker_log["sigma(x,x')i.....([m,rad])"] = (sigma_x,sigma_xp)
    tracker_log["sigma(y,y')i.....([m,rad])"] = (sigma_y,sigma_yp)
    tracker_log["sigma(Dphi,w)i....([rad,])"] = (sigma_Dphi,sigma_w)
    tracker_log["sigma(z,Dp2p)i......([m,])"] = (sigma_z,sigma_Dp2p)
    tracker_log['betax_i......[m]']           = betax_i
    tracker_log['betay_i......[m]']           = betay_i
    tracker_log['betaw_i....[rad]']           = betaw
    tracker_log['betaz_i..[m/rad]']           = betaz
    tracker_log['emitx_i........[m]']         = emitx_i
    tracker_log['emity_i........[m]']         = emity_i
    tracker_log['emitw_i,wmx..[rad]']         = (emitw, wmx)
    tracker_log['emitz_i........[m]']         = emitz
    tracker_log['Dp2p,Dp2pmx..[%]']           = (Dp2p0*1.e2,Dp2pmx*1.e2)
    tracker_log['acceptance Dp2p...[%]']      = PARAMS['Dp2pAcceptance']*1.e2
    tracker_log['accpetance z.....[mm]']      = PARAMS['zAcceptance']*1.e3
    tracker_log['lattice version......']      = PARAMS['lattice_version']
    tracker_log['mapping..............']      = PARAMS['mapping']
    tracker_log['DT/T-kin.............']      = PARAMS['DT2T']
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
    t2 = time.process_time()
    track_soll(lattice)  # <----- track soll
    t3 = time.process_time()
    # TmStamp.stamp('START TRACK')
    progress(('(track design)', '(track bunch)', '', ''))
    live_lost = track(lattice,bunch,options) # <----- track bunch returns (live,lost)-bunch
    t4 = time.process_time()
    # print(TmStamp.as_str())

    # make 2D projections
    if show:
        print('FILL PLOTS')
        projections(live_lost)
    t5 = time.process_time()
    if save:
        print('SAVE FRAMES')
        frames(lattice, skip)
    t6 = time.process_time()
    if losses:
        print('SHOW LOSSES')
        loss_plot(lattice,live_lost)

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

#----------------main------------
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
        input_file = 'yml/trackINstat.yml'
        if len(sys.argv) == 2:
            input_file    = sys.argv[1]
        print('input="{}"'.format(input_file))
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
    options['particles_per_bunch'] = 1000*1
    options['show']                = True
    options['save']                = False
    options['skip']                = 1
    options['losses']              = False

    # start the run
    tracker(options)
