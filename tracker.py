#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
___version___='v10.2.0'
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
#TODO: statistical analalysis of bunch: position, size, etc ...
#TODO: how to get the hokey stick?
#TODO: check w-acceptance at each node entrance
#TODO: no phase damping - why?
import sys,os
import numpy as NP
import matplotlib.pyplot as plt
import time
from string import Template
from math import sqrt, degrees, radians, ceil
import pprint, inspect
import argparse

from lattice_generator import factory
import elements as ELM
from setutil import PARAMS, FLAGS, dictprnt, Ktp, waccept
from setutil import WConverter, Functions
from bunch import BunchFactory, Gauss1D, Track, Tpoint, Bunch
import PoincareMarkerAgent as pcmkr

# max track amplitudes are set here!
xlim_max  = ylim_max  =  10.e-3
xplim_max = yplim_max =  10.e-3
zlim_max  = zplim_max = 100.e-3
limit     = \
    sqrt(xlim_max**2+xplim_max**2+ylim_max**2+yplim_max**2+zlim_max**2+zplim_max**2)

def PRINT_PRETTY(obj):
    file = inspect.stack()[0].filename
    print('DEBUG_ON ==============>  '+file)
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

def scatterPlot(live_lost, abszisse, ordinate, text, minmax=(1.,1.)):
    """ 
    The scatter plots 
    """
    live_bunch, lost_bunch = live_lost
    txt = ('IN {}'.format(text),'OUT {}'.format(text))
    initial = 0; final = -1

    # ENTRANCE-DATA
    loc = initial
    x=[]; y=[]; xlost=[]; ylost=[]
    nbprt = live_bunch.nbparticles()+lost_bunch.nbparticles()
    for particle in iter(live_bunch): # particles
        track  = particle.track
        tpoint = track.getpoints()[loc]
        point  = tpoint()
        x.append(point[abszisse]*1e3)
        y.append(point[ordinate]*1e3)
    xmax0 = max([abs(i) for i in x])
    ymax0 = max([abs(i) for i in y])
    
    if lost_bunch.nbparticles() != 0:     # lost particles
        for particle in iter(lost_bunch):
            track  = particle.track
            tpoint = track.getpoints()[loc]
            point  = tpoint()
            xlost.append(point[abszisse]*1e3)
            ylost.append(point[ordinate]*1e3)
        xmax1 = max([abs(i) for i in xlost])
        ymax1 = max([abs(i) for i in ylost])
    else:
        xmax1=xmax0
        ymax1=ymax0

    # Axis scales
    xmax = max(xmax0,xmax1)*1.1
    ymax = max(ymax0,ymax1)*1.1
    
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

    # EXIT-DATA
    loc = final
    x=[]; y=[]
    nbprt = live_bunch.nbparticles()
    for particle in iter(live_bunch): # particles
        track  = particle.track
        tpoint = track.getpoints()[loc]
        point  = tpoint()
        x.append(point[abszisse]*1e3)
        y.append(point[ordinate]*1e3)
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
    """ 2D phase space projection at marker position """
    plt.figure()    # new figure instance
    agent_cnt = 0
    frames = []
    # gather and count markers
    for node in iter(lattice):
        if isinstance(node,ELM.MRK) and isinstance(node.agent,pcmkr.PoincareMarkerAgent):
            marker_position = node.position
            marker_agent    = node.agent
            agent_cnt += 1
            if agent_cnt%skip == 0:
                frames.append((agent_cnt,marker_agent))
    # make an estimate for x- and y-axis
    dummy,first_frame = frames[0]
    x = [abs(tp()[first_frame.xaxis]) for tp in first_frame.tpoints]
    y = [abs(tp()[first_frame.yaxis]) for tp in first_frame.tpoints]
    xmax = max(x)*1.5
    ymax = max(y)*1.5
    # invoke actions on Marker
    for agent_cnt,agent in iter(frames):
        marker = agent.parent     # agent's parent is MRK object
        marker_position = marker.position
        marker.do_actions(agent_cnt,xmax,ymax,marker_position)
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()
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
    """ Show progress """
    template = Template('$tx1 $tx2 $tx3 $tx4')
    res = template.substitute(tx1=tx[0] , tx2=tx[1] , tx3=tx[2] , tx4=tx[3] )
    # print('\r{}\r'.format(res),end="")
    sys.stdout.write('{}\r'.format(res))
def track_node(node,particle,options):
    """ Tracks a particle through a node """
    def norm(track):
        track_norm = \
        sqrt(track[0]**2+track[1]**2+track[2]**2+track[3]**2+track[4]**2+track[5]**2)
        return track_norm > limit
    
    track   = particle.track
    last_tp = track.getpoints()[-1]
    lost    = False
    try:
        """ maping happens here! """
        new_point = node.map(last_tp())
        new_tp    = Tpoint(point=new_point)
    except (ValueError,OverflowError,ELM.OutOfRadialBoundEx) as ex:
        reason = ex.__class__.__name__
        s = last_tp()[8]
        print('@map in track_node: {} at s={:6.2f} [m]'.format(reason,s))
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
        if abs(new_point[Ktp.zp]) < PARAMS['Dp2p_Acceptance']:
            lost = False
        else:
            lost = True
        # check z-acceptance
        if abs(new_point[Ktp.z]) < PARAMS['z_Acceptance']:
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
        if isinstance(node,ELM.MRK) and isinstance(node.agent,pcmkr.PoincareMarkerAgent):
            node.agent.add_track_point(new_tp)
    particle.lost = lost
    return lost
def track(lattice,bunch,options):
    """
    Tracks a bunch of particles through the lattice using maps
    - lattice is a list of elements (class Node)
    - bunch (class Bunch) is a list of particles (class Particle)
    - each particle in a bunch has a track=particle['track']
    - track (class Track) is a list of points (class Tpoint)
    
    Input: lattice , bunch, options
    """
    ndcnt  = 0
    lnode  = len(lattice.seq)
    pgceil =  ceil(lnode/100)    # every 1% progress update
    nlost  = 0
    nbpart = bunch.nbparticles()
    lbunch = Bunch()    # lost particles go into this bunch
    printProgressBar(0,lnode,prefix="Progress:",suffix="Complete",length=50)
    for node in iter(lattice):              # nodes
        ndcnt +=1
        for particle in iter(bunch):        # particles
            lost = track_node(node,particle,options)
            if lost:
                lbunch.addparticle(particle)
                bunch.removeparticle(particle)
                nlost += 1
        # showing track-loop progress
        if ndcnt%pgceil == 0 or ndcnt == lnode: 
            printProgressBar(ndcnt,lnode,prefix="Progress:",suffix="Complete",length=50)
    live = nbpart - lbunch.nbparticles()
    print('\nTRACKING DONE (live particles {}, lost particles {})               '.format(live,nlost))
    return (bunch,lbunch)
def track_soll(lattice, injection_energy, start_position=0.):
    # TODO kann in der neuen Version wegfallen?
    """
    Track the reference particle through the lattice.
    NEW: Energy adjustemenr is already done when the node is added to the lattice.
    OLD: and redefines the lattice 
    element parameters according to the energy of the accelerated reference particle.
    """
    soll_track = Track()
    tp0 = Tpoint(NP.array([ 0., 0., 0., 0., 0., 0., injection_energy, 1., start_position, 1.]))
    soll_track.addpoint(tp0)   # 1st track point
    for node in iter(lattice):
        pi = soll_track.getpoints()[-1]   # track point at entrance
        """ OLD: energy adjustment """
        # node.adjust_energy(pi()[Ktp.T])
        # apply soll map 
        # pf = node.soll_map(pi())
        pf = node.map(pi())
        tpoint = Tpoint(pf)               # Tpoint at exit
        soll_track.addpoint(tpoint)
    return soll_track
def tracker(input_file,options):
    """  Prepare and launch tracking  """
    npart = options['particles_per_bunch']
    print('-----------------------track_bunch with {} particles---'.format(npart))

    # !!FIRST!! make lattice
    t0       = time.process_time()
    filepath = input_file
    lattice  = factory(filepath)

    # No tracking without acceleration
    if FLAGS['dWf'] == 0:
        print('{}'.format('IMPOSSIBLE: no tracking without acceleration!'))
        sys.exit()

    # calculate twiss paramters at entrance
    waccept(lattice.first_gap)
    tkin     = PARAMS['injection_energy']
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
    # Dp2pmx     = PARAMS['Dp2pmx']
    Dp2p0      = PARAMS['Dp2p0']
    # {Dphi,w}  T.Wangler units
    tww = PARAMS['twiss_w_i']
    betaw,alfaw,gammaw,emitw = tww()
    sigma_Dphi  = tww.sigmaH()
    sigma_w     = tww.sigmaV()
    wmax         = PARAMS['wmax']

    # gather for print
    tracker_log = {}
    tracker_log['T-kin..............[MeV]'] = tkin
    tracker_log["sigma(x,x')i...([m,rad])"] = (sigma_x,sigma_xp)
    tracker_log["sigma(y,y')i...([m,rad])"] = (sigma_y,sigma_yp)
    tracker_log["sigma(Dphi,w)i..([rad,])"] = (sigma_Dphi,sigma_w)
    tracker_log["sigma(z,Dp2p)i....([m,])"] = (sigma_z,sigma_Dp2p)
    tracker_log['betax_i..............[m]'] = betax_i
    tracker_log['betay_i..............[m]'] = betay_i
    tracker_log['betaw_i............[rad]'] = betaw
    tracker_log['betaz_i..........[m/rad]'] = betaz
    tracker_log['emitx_i..............[m]'] = emitx_i
    tracker_log['emity_i..............[m]'] = emity_i
    tracker_log['emitw_i,wmax.......[rad]'] = (emitw, wmax)
    tracker_log['emitz_i..............[m]'] = emitz
    # tracker_log['Dp2p,Dp2pmx..........[%]'] = (Dp2p0*1.e2,Dp2pmx*1.e2)
    tracker_log['Dp2p.................[%]'] = Dp2p0*1.e2
    tracker_log['acceptance Dp2p......[%]'] = PARAMS['Dp2p_Acceptance']*1.e2
    tracker_log['accpetance z........[mm]'] = PARAMS['z_Acceptance']*1.e3
    tracker_log['lattice version.........'] = PARAMS['lattice_version']
    tracker_log['mapping.................'] = PARAMS['mapping']
    tracker_log['DT/T-kin................'] = PARAMS['DT2T']
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
    track_soll(lattice, PARAMS['injection_energy'])  # <----- track soll
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
    track_soll(lattice, PARAMS['injection_energy'])  # <----- track soll
    table = sollTrack.as_table()
    DEBUG_TEST0('sollTrack:\n'+table)
    first = sollTrack[0]
    last  = sollTrack[-1]
    DEBUG_TEST0('sollTrack:\n(first): {}\n (last): {}'.format(first.as_str(),last.as_str()))
#----------------main------------
if __name__ == '__main__':
    DEBUG_TEST0 = DEBUG_OFF
    # use ArgumentParser to put result in 'args'
    parser = argparse.ArgumentParser(prog='python tracker.py')
    group  = parser.add_mutually_exclusive_group()
    group1 = parser.add_mutually_exclusive_group()
    parser.add_argument("--p", metavar="N", default=1750, type=int,   help="N particles per bunch")
    parser.add_argument("--hide", action="store_false",               help="hide IN/OUT scatter plots")
    group.add_argument ("--file", default="trackerINwork.yml",        help="lattice input-file")
    group.add_argument ("--tmpl",                                     help="template number")
    parser.add_argument("--run",                                      help="run number")
    group1.add_argument("--pcuts", action="store_true",               help="save poincare cuts")
    group1.add_argument("--losses", action="store_true",              help="run in losses mode")
    parser.add_argument("--skip", metavar="N", default="1", type=int, help="skip every N poincare cuts")
    args = vars(parser.parse_args())
    # DEBUG_ON(args)
    options = {}
    options['particles_per_bunch'] = args['p']
    options['show']                = args['hide']
    options['save']                = args['pcuts']
    options['skip']                = args['skip']+1
    options['losses']              = args['losses']

    print('tracker.py {} on python {}.{}.{} on {}'.format(___version___,sys.version_info.major,sys.version_info.minor,sys.version_info.micro,sys.platform))
    
    # adapt to legacy code which uses 'Args'
    Args  = {}
    tmpl  = args['tmpl']
    run   = args['run'] 
    Args['mode']  = 'no_m4' if tmpl == None else 'm4'
    Args['file']  = args['file']
    Args['tmpl']  = ''
    Args['macro'] = ''
    if Args['mode'] == 'm4':
        Args['tmpl']   = 'yml/tmpl_{}.yml'.format(tmpl)
        Args['macro']  = 'yml/macros_{}.{}.sh'.format(tmpl,run) if run != None else 'yml/macros_{}.sh'.format(tmpl)
    print('This run: input({}), template({}), macro({})'.format(Args['file'],Args['tmpl'],Args['macro']))

    # let's go. All  input is parsed...
    input_file = Args['file']
    if sys.platform == 'win32':
        if Args['mode']   == 'no_m4':
            pass
        elif Args['mode'] == 'm4':
            command = 'yml\m4_launch.bat {} {} {}'.format(Args['file'],Args['tmpl'],Args['macro'])
            stat = os.system(command)
            if stat != 0:
                print('\nWARNING: system-command returned error - using standard "yml/simuIN.yml" without m4-preprocessing!')
                print(  'WARNING: system-command returned error - using standard "yml/simuIN.yml" without m4-preprocessing!')
                print(  'WARNING: system-command returned error - using standard "yml/simuIN.yml" without m4-preprocessing!\n')
        else:
            print('Internal error!')
            sys.exit(1)
    elif sys.platform == 'darwin' or sys.platform.startswith('linux'):
        if Args['mode']   == 'no_m4':
            pass
        elif Args['mode'] == 'm4':
            macros_file   = Args['macro']
            template_file = Args['tmpl']
            # launch macros script with bash
            command = 'chmod +x {}'.format(macros_file)
            command = "{0};{1} {2} {3}".format(command,macros_file,template_file, input_file)            
            os.system(command)
        else:
            print('Internal error!')
            sys.exit(1)
    else:
        print('wrong platform')
        sys.exit(1)

    # start the tracking
    tracker(input_file,options)
