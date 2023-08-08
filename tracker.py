#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
__version__='v11.0.2'
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
#TODO: how to get the hokey stick? - done
#TODO: check w-acceptance at each node entrance - done
#TODO: no phase damping - why? - solved with version 10.0.0

import sys,os
import numpy as np
import matplotlib.pyplot as plt
import time
from math import sqrt, degrees, radians, ceil,pi
import argparse
import unittest
import h5py


from lattice_generator import factory
import elements as ELM
from setutil import PARAMS, FLAGS, dictprnt, Ktp, WConverter
from setutil import RUN_MODE, Functions, DEBUG_ON, DEBUG_OFF,Proton
from bunch import BunchFactory, Gauss1D, Gauss2D, Track, Tpoint, Bunch
from PoincareMarkerAgent import PoincareMarkerAgent
from trackPlot import scatter11,scatterInOut

# max limits for track amplitudes
xlim_max  = 0.1
ylim_max  = 0.1
zlim_max  = 0.1
limit_xyz = sqrt(xlim_max**2+ylim_max**2+zlim_max**2)

class Fifo:
    def __init__(self):
        self._first = None
        self._last = None
        self._max = 10000
        self._cnt = 0
    def append(self,data):
        self._cnt += 1
        if self._cnt > self._max: 
            return
        node = [data,None]
        if self._first is None:
            self._first = node
        else:
            self._last[1] = node
        self._last = node
    def pop(self):
        if self._first is None:
            return None
        node = self._first
        self._first = node[1]
        return node[0]
    @property
    def max(self):
        return self._max
    @max.setter
    def max(self,value):
        self._max = value
# txt FIFOs
fifo    = Fifo()   # halo-losses
fifo_m  = Fifo()   # map-losses
fifo_z  = Fifo()   # z-losses
fifo_xy = Fifo()   # xy-losses
# position FIFOs
sfifo    = Fifo()
sfifo_xy = Fifo()
sfifo_m  = Fifo()
sfifo_z  = Fifo()

def projections_1(lattice,live_lost):
    """ 2D projections of 6D phase space """
    
    def projection(live_lost,ix,iy,fig_txt,scale=[(1.,1.),(1.,1.)]):
        """ projection: scatter plots at IN and OUT 
            ix: track coordinate of plot abscissa
            iy: track coordinate of plot ordinate
            fig_text: annotation
            live_lost: tuple of particle containers with (live,lost)-particles in bunch
            scale: scales for IN and OUT plots [(xIN,yIN),(xOUT,yOUT)]
        """
        inout=dict(IN=0,OUT=-1)
        live,lost=live_lost
        nblive=live.nbparticles()
        nblost=lost.nbparticles()

        golden = (1.+sqrt(5.))/2.; width = 10; height = width/golden
        plt.figure(num=fig_txt,constrained_layout=False, figsize=(width, height))

        xlive=np.zeros(2)
        ylive=np.zeros(2)
        for particle in iter(live):
            track=particle.track
            point=track.getpoints()[inout['IN']]()
            xlive=np.append(xlive,point[ix]*scale[inout['IN']][0])
            ylive=np.append(ylive,point[iy]*scale[inout['IN']][1])
        livemax=np.array([np.amax(np.abs(xlive)),np.amax(np.abs(ylive))])
        xloss=np.zeros(2)
        yloss=np.zeros(2)
        for particle in iter(lost):
            track=particle.track
            point=track.getpoints()[inout['IN']]()
            xloss=np.append(xloss,point[ix]*scale[inout['IN']][0])
            yloss=np.append(yloss,point[iy]*scale[inout['IN']][1])
        lossmax=np.array([np.amax(np.abs(xloss)),np.amax(np.abs(yloss))])
        xymax=np.fmax(livemax,lossmax)*1.03 # add 3% margin
        box_txt=f'IN {fig_txt} {nblive+nblost} particles'

        xlive1=np.zeros(2)
        ylive1=np.zeros(2)
        for particle in iter(live):
            track=particle.track
            point=track.getpoints()[inout['OUT']]()
            xlive1=np.append(xlive1,point[ix]*scale[inout['OUT']][0])
            ylive1=np.append(ylive1,point[iy]*scale[inout['OUT']][1])
        livemax1=np.array([np.amax(np.abs(xlive1)),np.amax(np.abs(ylive1))])
        xloss1=np.zeros(2)
        yloss1=np.zeros(2)
        box1_txt=f'OUT {fig_txt} {nblost} lost particles'

        plotmax=np.array([max(xymax[0],livemax1[0]),max(xymax[1],livemax1[1])])
        DEBUG_OFF(f'plotmax={plotmax}')
        ax=plt.subplot(121)
        scatterInOut(xlive,ylive,xloss,yloss,plotmax,box_txt,ax)
        ax=plt.subplot(122)
        scatterInOut(xlive1,ylive1,xloss1,yloss1,plotmax,box1_txt,ax)
        return 

    # TODO can this be made more elegantly? This function made me much pain.
    def projection_dPhidW(live_lost,lattice):
        """
        To get scales for {Dphi-DW} plots from internal {Dz-Dp/p} coordinates
        we need the reference energy of each particle and the gap's frequency.
        We take the gap's frequency from the cavity.
        We take the reference energy from a bunch particle that survived.
        """
        DELTA='\u0394'
        PHI  ='\u03A6'
        GAMMA='\u03b3'

        live,lost=live_lost
        # first,last gap
        in_gap= lattice.first_gap
        out_gap=lattice.last_gap
        DEBUG_OFF(f'1st gap: {in_gap.toString()}')
        DEBUG_OFF(f'last gap: {out_gap.toString()}')
        # frequencies of first,last
        freqIN= in_gap.freq
        freqOUT=out_gap.freq
        # trak-points of particle 0 in bunch
        points=live.getparticles()[0].track.getpoints()
        # first track-point
        point=points[0]
        # kin energy of first track point
        tkIN=point()[Ktp.T]
        # last track-point
        point=points[-1]
        # kin energy of last track-point
        tkOUT=point()[Ktp.T]
        DEBUG_OFF(f'freq(IN,OUT) {(freqIN,freqOUT)}  tk(IN,OUT) {(tkIN,tkOUT)}')
        DEBUG_ON(f'(W-IN,W-OUT)={(tkIN,tkOUT)}')
        convIN =WConverter(tkIN,freqIN)
        convOUT=WConverter(tkOUT,freqOUT)
        # scale=[
        #      (degrees(convIN.zToDphi(1.)),   convIN.Dp2pTow(1.)*1e2),
        #      (degrees(convOUT.zToDphi(1.)),  convOUT.Dp2pTow(1.)*1e2)
        #       ]
        # scale=[
        #      (degrees(convIN.zToDphi(1.)),   convIN.Dp2pToDW2W(1.)*1e2),
        #      (degrees(convOUT.zToDphi(1.)),  convOUT.Dp2pToDW2W(1.)*1e2)
        #       ]
        scale=[
             (degrees(convIN.zToDphi(1.)),   convIN.Dp2pToDW(1.)*1e3),
             (degrees(convOUT.zToDphi(1.)),  convOUT.Dp2pToDW(1.)*1e3)
              ]
        return projection(live_lost,Ktp.z,Ktp.zp,f'{DELTA}{PHI}-{DELTA}W [deg,KeV]',scale=scale)
    
    DELTA='\u0394'
    # longitudinal
    # projection(live_lost,Ktp.z,Ktp.zp,f'z-{DELTA}p/p [m,]')
    # projection(live_lost,Ktp.z,Ktp.zp,f'z-{DELTA}p/p [mm,%]',scale=[(1.e3,1.e2),(1.e3,1.e2)])
    projection_dPhidW(live_lost,lattice)
    # transverse
    # projection(live_lost,Ktp.x,Ktp.y, f"x-y [mm,mrad]",scale=[(1.e3,1.e3),(1.e3,1.e3)])
    projection(live_lost,Ktp.x,Ktp.xp,f"x-x' [mm,mrad]",scale=[(1.e3,1.e3),(1.e3,1.e3)])
    projection(live_lost,Ktp.y,Ktp.yp,f"y-y' [mm,mrad]",scale=[(1.e3,1.e3),(1.e3,1.e3)])
    plt.show()
    return
def frames(lattice, skip):
    """ 2D phase space projection at marker position """
    plt.figure()    # new figure instance
    agent_cnt = 0
    frames = []
    # gather and count markers

    for node in iter(lattice):
        if isinstance(node,PoincareMarkerAgent):
            agent = node
            agent_cnt += 1
            if agent_cnt%skip == 0:
                frames.append((agent_cnt,agent))
    # make an estimate for x- and y-axis
    lrx = options['lrx']
    dummy,lrx_frame = frames[lrx]
    x = [abs(tp()[lrx_frame.xaxis]) for tp in lrx_frame.tpoints]
    y = [abs(tp()[lrx_frame.yaxis]) for tp in lrx_frame.tpoints]
    xmax = max(x)*1.5
    ymax = max(y)*1.5
    # invoke actions on Marker
    for agent_cnt,agent in iter(frames):
        position = agent.position
        agent.do_action(agent_cnt,xmax,ymax,position)
def progress_bar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
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
def loss_histograms(lattice,fifos,binsize=5):
    """ make histogram plots of losses captured in FIFO buffers """
    latlen = lattice.length
    bins   = int(latlen/binsize)
    golden = (1.+sqrt(5.))/2.; width = 10; height = width/golden
    fig    = plt.figure(num='losses',figsize=(width, height),layout='constrained')
    axs    = fig.subplots(2,2)
    for ff_cnt,fifo in enumerate(fifos):
        fig_txt   = fifo['title']
        fifo_data = fifo['fifo']
        sdata = []
        while(True):
            data = fifo_data.pop()
            if data is None: break
            sdata.append(data)
        ax = axs.flat[ff_cnt]
        ax.set_title(fig_txt)
        ax.set_xlabel(r"s[m]")
        ax.set_ylabel(r"# particles")
        ax.hist(sdata,bins,range=(0.,latlen))
    plt.show()
def track_node_1(node,particle,options):
    """ Tracks a particle through a node """
    def norm(tp):
        tpnorm = sqrt(tp()[Ktp.x]**2+tp()[Ktp.y]**2+tp()[Ktp.z]**2)
        return tpnorm > limit_xyz
    
    track   = particle.track
    last_tp = track.getpoints()[-1]
    s       = last_tp()[Ktp.S]
    lost    = False

    # ********************************************************************************
    # maping happens here! (map-losses)
    try:
        new_point = node.map(last_tp())
        new_tp    = Tpoint(point=new_point)
    except (ValueError,OverflowError,ELM.OutOfRadialBoundEx) as ex:
        txt = ex.message
        DEBUG_OFF(txt)
        fifo_m.append(txt)
        sfifo_m.append(s)
        lost = True
        particle.lost = lost
        return lost
     # ********************************************************************************
   
    # check new_tp against reasonable physical limits (halo-losses)
    if norm(new_tp):
        fifo.append(f'halo limits at {s:.4e} m')
        sfifo.append(s)
        lost = True
        particle.lost = lost
        return lost
    
    # aperture checks (z-losses and xy-losses)
    if FLAGS['useaper']:
        lost = node.aper_check(new_tp,s,fifo_z=fifo_z,sfifo_z=sfifo_z,fifo_xy=fifo_xy,sfifo_xy=sfifo_xy)
    particle.lost = lost
    # eins=1
    # replace old tp by new tp (save memory)
    if track.nbpoints() > 1:
        track.removepoint(last_tp)
    track.addpoint(new_tp)
    # PoincareMarker keeeps all tps to dump frames
    if isinstance(node,PoincareMarkerAgent):
        node.appendPhaseSpace(new_tp)
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
    if options['h5dump']:
        h5File_name = options['h5file']
        try:
            os.remove(h5File_name)
        except OSError as e:
            pass
        finally:
            print(f'creating new HDF5-file "{h5File_name}"')
            h5File=h5py.File(h5File_name,"w-")
        h5frames_grp=h5File.create_group('/frames')

    nb_nodes     = len(lattice.seq)
    nb_particles = bunch.nbparticles()
    pgceil       =  ceil(nb_nodes/100)    # every 1% progress update
    lbunch       = Bunch()    # lost particles go into this bunch
    progress_bar(0,nb_nodes,prefix="Progress:",suffix="Complete",length=50)

    for n_cnt,node in enumerate(iter(lattice)):              # nodes
        # current number of particles in bunch
        current_nb_particles = bunch.nbparticles()
        h5dump = options['h5dump'] and (n_cnt%options['h5skip'] == 0)
        # HDF5 dumping: create data set
        if h5dump: 
            h5ds = h5frames_grp.create_dataset(f'{n_cnt}',(current_nb_particles,10),dtype='f8')
        for p_cnt,particle in enumerate(iter(bunch)):        # particles
            lost = track_node_1(node,particle,options)
            if lost:
                lbunch.addparticle(particle)
                bunch.removeparticle(particle)
            else:
                # HDF5 dumping: fill data set
                if h5dump:
                    tp = particle.track.getpoints()[-1]()
                    DEBUG_OFF(f'node={n_cnt} particle={p_cnt} track-point={tp}')
                    h5ds[p_cnt] = tp
                pass
        # showing track-loop progress
        if n_cnt%pgceil == 0 or n_cnt == nb_nodes: 
            progress_bar(n_cnt,nb_nodes,prefix="Progress:",suffix="complete",length=50)

    DEBUG_OFF(f'n_cnt= {n_cnt}')
    lost = lbunch.nbparticles()
    print('\nTRACKING DONE (particles {}, live {}, lost {})'.format(nb_particles,nb_particles-lost,lost))
    # closing HDF5 file
    if options['h5dump']:
        h5File.close()
    return (bunch,lbunch)
def tracker(input_file,options):
    """ Prepare and launch tracking  """
    # fifo limits
    fifo_m.max  = 10
    fifo_z.max  = 10
    fifo_xy.max = 10

    npart = options['particles_per_bunch']
    # make lattice
    t0       = time.process_time()
    filepath = input_file
    lattice  = factory(filepath)
    DEBUG_OFF(PARAMS['twiss_w_i']())
    # w acceptance
    FLAGS['accON'] = lattice.accON
    DEBUG_OFF(PARAMS['twiss_w_i']())
    if not FLAGS['accON']:
        # no acceleration
        print('{}'.format('IMPOSSIBLE: no tracking without acceleration!'))
        sys.exit()
    # run_mode
    FLAGS['mode'] = RUN_MODE[1]
    print(f'running in \'{FLAGS["mode"]}\' mode')

    t1 = time.process_time()
    # bunch-configuration from PARAMS
    # {x,xp}   standard units
    twx = PARAMS['twiss_x_i']
    betax_i,alfax_i,gammax_i,emitx_i = twx()
    sigma_x   = twx.sigmaH()
    sigma_xp  = twx.sigmaV()
    DEBUG_OFF(f'{{x}}x{{xp}} {twx()}')

    # {y,yp}
    twy = PARAMS['twiss_y_i']
    betay_i,alfay_i,gammay_i,emity_i = twy()
    sigma_y   = twy.sigmaH()
    sigma_yp  = twy.sigmaV()
    DEBUG_OFF(f'{{y}}x{{yp}} {twy()}')

    # {z,Dp2p}  T3D units
    twz = PARAMS['twiss_z_i']
    betaz_i,alfaz_i,gammaz_i,emitz_i = twz()
    sigma_z    = twz.sigmaH()
    sigma_Dp2p = twz.sigmaV()
    DEBUG_OFF(f'{{z}}x{{Dp2p}} {twz()}')
    Dp2p0      = PARAMS['Dp2p0_i']

    # {Dphi,w}  T.Wangler units
    tww = PARAMS['twiss_w_i']
    betaw_i,alfaw_i,gammaw,emitw_i = tww()
    sigma_Dphi  = tww.sigmaH()
    sigma_w     = tww.sigmaV()
    DEBUG_OFF(f'{{Dphi}}x{{w}} {tww()}')

    # gather for print
    tracker_log = {}
    tracker_log['Description.............']           = PARAMS.get('descriptor')
    tracker_log['mapping.................']           = FLAGS['mapping']
    tracker_log['useaper.................']           = FLAGS['useaper']
    tracker_log['Options.................']           = f'{options}'
    tracker_log['Tk_i...............[MeV]']           = '{} kin. energy @ injection'.format(lattice.injection_energy)
    tracker_log['acceptance..\u0394p/p.....[%]']      = f"{PARAMS['Dp2pmax']*1.e2:.3f}"
    tracker_log['acceptance..\u0394\u03B3..........'] = f"{PARAMS['wmax']:.2e}"
    tracker_log['\u03B2w_i...............[rad]']      = betaw_i
    tracker_log['\u03B2x_i.................[m]']      = betax_i
    tracker_log['\u03B2y_i.................[m]']      = betay_i
    tracker_log['\u03B2z_i.............[m/rad]']      = betaz_i
    tracker_log['\u03B5w_i..{\u0394\u03C6,\u0394\u03B3}......[rad]'] = f"{emitw_i:.2e}"
    tracker_log['\u03B5x_i.................[m]']      = emitx_i
    tracker_log['\u03B5y_i.................[m]']      = emity_i
    tracker_log['\u03B5z_i..{z,\u0394p/p}.......[m]'] = emitz_i
    tracker_log['lattice version.........']           = PARAMS['lattice_version']
    tracker_log["\u03C3(x,x')_i......([m,rad])"]      = "{:.2e} {:.2e}".format(sigma_x,sigma_xp)
    tracker_log["\u03C3(y,y')_i......([m,rad])"]      = "{:.2e} {:.2e}".format(sigma_y,sigma_yp)
    tracker_log["\u03C3(z,\u0394p/p)_i.......([m,])"] = "{:.2e} {:.2e}".format(sigma_z,sigma_Dp2p)
    tracker_log["\u03C3(\u0394\u03C6,\u0394\u03B3)_i......([rad,])"] = "{:.2e} {:.2e}".format(sigma_Dphi,sigma_w)
    tracker_log["\u03C3(\u0394\u03C6,\u0394\u03B3)_i......([deg,])"] = "{:.2e} {:.2e}".format(degrees(sigma_Dphi),sigma_w)
    tracker_log['\u0394p/p0................[%]']      = f"{Dp2p0*1.e2:.3f}"
    tracker_log['\u0394T/T_i..(\u0394W/W).......[%]'] = f"{PARAMS['DT2T_i']*1e2:.3f} kin. energy spread"
    dictprnt(tracker_log,'Tracker Log',njust=36); print()

    # bunch factory
    bunchfactory = BunchFactory()
    bunchfactory.twiss=(twx,twy,twz)  # twiss parameters for initial bunch distribution
    bunchfactory.numberOfParticles=npart
    bunch = bunchfactory()
    t2 = t3 = time.process_time()

    # ********************************************************************************
    # track: returns tuple (live,lost) bunches
    live_lost = track(lattice,bunch,options) 
    t4 = time.process_time()
    # ********************************************************************************

    # make 2D projections
    if options['show']:
        print('FILL PLOTS')
        projections_1(lattice,live_lost)
    t5 = time.process_time()
    
    if options['save']:
        print('SAVE FRAMES')
        frames(lattice, options['skip'])
    t6 = time.process_time()
    
    if options['losses']:
        print('SHOW LOSSES')
        fifos = [
            dict(title='transverse losses',fifo=sfifo_xy),
            dict(title='longitudinal losses',fifo=sfifo_z),
            dict(title='mapping losses',fifo=sfifo_m),
            dict(title='halo losses',fifo=sfifo),
            ]   
        loss_histograms(lattice,fifos)
    t7 = time.process_time()

    # finish up
    print()
    ttotal = t7-t0
    print('total time     >> {:6.3f} [sec]'.format(ttotal))
    print('parse lattice  >> {:6.3f} [sec] {:4.1f} [%]'.format((t1-t0),(t1-t0)/(ttotal)*1.e2))
    print('generate bunch >> {:6.3f} [sec] {:4.1f} [%]'.format((t2-t1),(t2-t1)/(ttotal)*1.e2))
    print('track bunch    >> {:6.3f} [sec] {:4.1f} [%]'.format((t4-t3),(t4-t3)/(ttotal)*1.e2))
    print('fill plots     >> {:6.3f} [sec] {:4.1f} [%]'.format((t5-t4),(t5-t4)/(ttotal)*1.e2))
    print('save frames    >> {:6.3f} [sec] {:4.1f} [%]'.format((t6-t5),(t6-t5)/(ttotal)*1.e2))
    print('bin losses     >> {:6.3f} [sec] {:4.1f} [%]'.format((t7-t6),(t7-t6)/(ttotal)*1.e2))

    while True:
        data = fifo.pop()
        if data is None: break
        DEBUG_OFF(data)
    while True:
        data = fifo_xy.pop()
        if data is None: break
        DEBUG_OFF(data)
    while True:
        data = fifo_m.pop()
        if data is None: break
        DEBUG_OFF(data)
    while True:
        data = fifo_z.pop()
        if data is None: break
        DEBUG_OFF(data)

class TestTracker(unittest.TestCase):
    def test_tracking(self):
        print('---------------------test_tracking---')
        print('what?')

#----------------main------------
if __name__ == '__main__':
    DEBUG_OFF(sys.argv)
    # ArgumentParser puts result in 'args'
    parser = argparse.ArgumentParser()
    group  = parser.add_mutually_exclusive_group()
    group1 = parser.add_mutually_exclusive_group()
    parser.add_argument("--p",      metavar="N", default=1750, type=int,     help="N particles per bunch")
    parser.add_argument("--hide",   action="store_true",                     help="hide IN/OUT scatter plots")
    group.add_argument ("--file",   default="trackerIN.yml",                 help="lattice input-file")
    group1.add_argument("--losses", action="store_true",                     help="run in losses mode")
    group1.add_argument("--pcuts",  action="store_true",                     help="save poincare cuts")
    parser.add_argument("--skip",   metavar="N", default="1", type=int,      help="skip every N poincare cuts")
    parser.add_argument("--lrx",    metavar="N", default="-1", type=int,     help="take N-th frame as axis limits. first=0, last=-1")
    parser.add_argument("--h5dump", action="store_true",                     help="dump tracks to HDF5 file")
    parser.add_argument("--h5skip", metavar="N", default="500", type=int,    help="skip every N track dumps")
    parser.add_argument("--h5file", default="frames.h5",                     help="HDF5 dump-file")
    args = vars(parser.parse_args())
    options = {}
    options['particles_per_bunch'] = args['p']
    options['show']                = not args['hide']
    options['save']                = args['pcuts']
    options['skip']                = args['skip']
    options['losses']              = args['losses']
    options['lrx']                 = args['lrx']
    options['h5dump']              = args['h5dump']
    options['h5skip']              = args['h5skip']
    options['h5file']              = args['h5file']

    # manipulate options
    if options['save']:
        options['show']   = False
        options['losses'] = False
    print('tracker.py {} on python {}.{}.{} on {}'.format(__version__,sys.version_info.major,sys.version_info.minor,sys.version_info.micro,sys.platform))
    
    # let's go. All  input is parsed...
    input_file = args['file']
    print('This run: input({})'.format(input_file))
    if False:
        # unittest.main()
        TestTracker(methodName='test_tracking').run()
    else:
        if sys.platform == 'darwin' or sys.platform.startswith('linux') or sys.platform == 'win32':
            print(f'Platform {sys.platform}')
        else:
            print('wrong platform')
            sys.exit(1)
        # start the tracking
        tracker(input_file,options)
