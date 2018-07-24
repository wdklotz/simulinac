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
import numpy as NP
import matplotlib.pyplot as plt
import time
from string import Template
from math import sqrt
# from joblib import Parallel, delayed

from lattice_generator import factory
import elements as ELM
from setutil import DEBUG, PARAMS, dictprnt, sigmas, K, PARAMS
from bunch import Track, Bunch, Gauss1D
from trackPlot import poincarePlot

# DEBUGGING
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass

def progress(tx):
    pass
    res = template.substitute(tx1=tx[0] , tx2=tx[1] , tx3=tx[2] , tx4=tx[3] )
    print('{}'.format(res), end="\r")

def scatterPlot(bunch, poincare_section, ordinate, abzisse, text, minmax):
    """ prepare the plot of a Poincaré section """
    nbsections = bunch.nbsections        # poincare sections per track
    psec = poincare_section
    if psec == 0:
        txt = '{} initial'.format(text)
    elif psec == (nbsections-1):
        txt ='{} final'.format(text)
    else:
        txt ='{} marker {}'.format(text, psec)
    fig = plt.figure(psec)
    x=[]; y=[]
    for t in bunch.tracks:                     # loop tracks
        x.append(t.point_at(psec)[ordinate]*1.e3)
        y.append(t.point_at(psec)[abzisse]*1.e3)
    good = (x,y)

    xi=[]; yi=[]
    for t in bunch.invalid_tracks:             # loop invalid tracks
        if psec < t.nbpoints:
            xi.append(t.point_at(psec)[ordinate]*1.e3)
            yi.append(t.point_at(psec)[abzisse]*1.e3)
    bad = (xi,yi)

    boxtext = '{} {} particles'.format(txt, bunch.nbtracks)
    poincarePlot(good, boxtext, minmax, bad = bad, projections = (1,1))
    return fig
    
def process_single_track(arg):
    """ successive mappings of a track through the lattice """
    ptrack  = arg[0]
    lattice = arg[1]
    invalid = False

    ti = ptrack.first()
    for element in lattice.seq:                  # loop elements
        try:
            tf = element.map(ti)                 # map!
        except ValueError as ex:
            invalid = True
            break
        ti = tf
        # take a poincare section at MRK elements only
        if isinstance(element,ELM.MRK) and not invalid:
            ptrack.append(tf)
    return (ptrack,invalid)

def track(lattice,bunch,smp=False):
    """
    Tracks a bunch of particles through the lattice using maps
    - lattice is a list of elements
    - bunch   is a list of particle tracks
    - smp flag False means do not use multipocessing
    """
    if DEBUG_TRACK == DEBUG_ON: dictprnt(bunch._params,'bunch',filter='tracklist')
    
    invalid_tracks = []
    valid_tracks   = []
    losses         = 0
    zeuge          = ('*\u007C*','**\u007C','*\u007C*','\u007C**')  # *|*
    tx4            = ' {}/{}/{} done/lost/initial'.format(0,0,bunch.nbtracks)

    if(smp):
        pass
        # arg = [(ptrack,lattice) for ptrack in bunch.tracks]
        # print()
        # r = Parallel(n_jobs=8, verbose=5)(map(delayed(process_single_track),arg))
        # trck,invalid = zip(*r)
        # for i,t in enumerate(trck):
        #     if invalid[i]:
        #         invalid_tracks.append(t)
        #     else:
        #         valid_tracks.append(t)
    else:
        for (tcount, ptrack) in enumerate(bunch.tracks):       # loop tracks
            invalid = process_single_track((ptrack, lattice))[1]
            if invalid:
                invalid_tracks.append(ptrack)
            else:
                valid_tracks.append(ptrack)
            # showing track-loop progress
            if (tcount+1)%25 == 0:
                losses = len(invalid_tracks)
                tx4    = ' {}/{}/{} done/lost/initial'.format(tcount+1, losses, bunch.nbtracks)
            progress(('(track design)', '(track bunch)', zeuge[tcount%4], tx4))
    # keep valid/invalid tracks in the bunch
    bunch.tracks         = valid_tracks
    bunch.invalid_tracks = invalid_tracks
    print('\ndone: good tracks {}, bad tracks {}'.format(bunch.nbtracks,  bunch.nbinvalid_tracks))
    return bunch

def track_soll(lattice):
    """
    Tracks the reference particle through the lattice and redefines the lattice element parameters to
    adapted to the energy of the accelerated reference particle.
    """
    soll_track = Track(start=NP.array([ 0., 0., 0., 0., 0., 0., PARAMS['sollteilchen'].tkin, 1., 0., 1.]))
    for element in lattice.seq:
        # DEBUG_SOLL_TRACK(element,' pos {:.4f} label "{}"'.format(element.position[1],element.label))
        ti = soll_track.last() #track: at entrance
        # DEBUG_SOLL_TRACK('track_soll(i) ', soll_track.point_str(ti))
        # DEBUG_SOLL_TRACK('track_soll: complete track\n{}'.format(soll_track.points_str()))
        """ energy adjustment """
        element.adjust_energy(ti[K.T])
        # DEBUG_SOLL_TRACK('track_soll: complete track\n{}'.format(soll_track.points_str()))
        """ mapping with soll map """
        tf = element.soll_map(ti)
        # DEBUG_SOLL_TRACK('track_soll(f) ',soll_track.point_str(tf))
        soll_track.append(tf)
        # DEBUG_SOLL_TRACK('track_soll: complete track\n{}'.format(soll_track.points_str()))
    
    # DEBUG_SOLL_TRACK('track_soll: complete track\n{}'.format(soll_track.points_str()))
    # DEBUG_SOLL_TRACK('track_soll(first) {}'.format(soll_track.first_str()))
    # DEBUG_SOLL_TRACK('track_soll(last)  {}'.format( soll_track.last_str()))
    return soll_track

def tracker(options):
    """ prepare and launch tracking """
    filepath          = options['filepath']
    particlesPerBunch = options['particlesPerBunch']
    show              = options['show']
    save              = options['save']
    skip              = options['skip']

    #make lattice with time
    t0 = time.clock()
    lattice = factory(filepath)
    t1 = time.clock()

    options['tkin [MeV'] = PARAMS['sollteilchen'].tkin
    dictprnt(options,'Tracker Options')
    print()

    # bunch-configuration
    alfax_i   = PARAMS['alfax_i'] 
    betax_i   = PARAMS['betax_i']
    emitx_i   = PARAMS['emitx_i']
    alfay_i   = PARAMS['alfay_i']
    betay_i   = PARAMS['betay_i']
    emity_i   = PARAMS['emity_i']
    sigmaz_i  = sqrt(PARAMS['emitz_i']*PARAMS['betaz_i'])   # Dz in [m]
    sigmazp_i = PARAMS['emitz_i']/sigmaz_i                  # Dp/p -n [-]

    sigmas_x  = sigmas(alfax_i, betax_i, emitx_i)
    sigmas_y  = sigmas(alfay_i, betay_i, emity_i)

    bunch = Bunch()
    bunch.nbtracks    = particlesPerBunch
    bunch.coord_mask  = (1,1,1,1,1,1)
    bunch.disttype    = Gauss1D
    bunch['sigma-x']  = sigmas_x[0]
    bunch['sigma-xp'] = sigmas_x[1]
    bunch['sigma-y']  = sigmas_y[0]
    bunch['sigma-yp'] = sigmas_y[1]
    bunch['sigma-z']  = 1.e-3
    bunch['sigma-zp'] = 1.e-3
    bunch['tkin']     = 70.
    bunch.populate_phase_space()

    # launch tracking and show final with time
    progress(('(track design)', '', '', ''))
    t2 = time.clock()
    track_soll(lattice)  # track soll
    t3 = time.clock()
    progress(('(track design)', '(track bunch)', '', ''))
    track(lattice,bunch) # track bunch
    t4 = time.clock()
    # make 2D projections
    # project_onto_plane(bunch, K.x, K.xp, show, save, skip)
    # project_onto_plane(bunch, K.x, K.y, show, save, skip)
    project_onto_plane(bunch, K.z, K.zp, show, save, skip)
    t5 = time.clock()

    print()
    print('total time     >> {:6.3f} [sec]'.format(t5-t0))
    print('parse lattice  >> {:6.3f} [sec] {:4.1f} [%]'.format((t1-t0),(t1-t0)/(t5-t0)*1.e2))
    print('generate bunch >> {:6.3f} [sec] {:4.1f} [%]'.format((t2-t1),(t2-t1)/(t5-t0)*1.e2))
    print('track design   >> {:6.3f} [sec] {:4.1f} [%]'.format((t3-t2),(t3-t2)/(t5-t0)*1.e2))
    print('track bunch    >> {:6.3f} [sec] {:4.1f} [%]'.format((t4-t3),(t4-t3)/(t5-t0)*1.e2))
    print('fill plots     >> {:6.3f} [sec] {:4.1f} [%]'.format((t5-t4),(t5-t4)/(t5-t0)*1.e2))

def project_onto_plane(bunch, ordinate, abzisse, show, save, skip):
    """ 2D phase space projections of Poincaré sections """
    symbol    = ("x","x'","y","y'","z","z'")
    sigmas    = (bunch['sigma-x'],bunch['sigma-xp'],bunch['sigma-y'],bunch['sigma-yp'],bunch['sigma-z'],bunch['sigma-zp'])
    text      = '{}-{}'.format(symbol[ordinate],symbol[abzisse])
    minmax    = (10.e3*sigmas[ordinate],10.e3*sigmas[abzisse])
    sections  = [i for i in range(bunch.nbsections)]
    for i in sections:
        if i != 0 and i != sections[-1] and i%skip != 0: continue   # skip sections, keep first and last
        fig = scatterPlot(bunch=bunch, poincare_section=i, ordinate=ordinate, abzisse=abzisse, text=text, minmax=minmax)
        if save:
            plt.savefig('figures/poincare_section_{}_{}.png'.format(text,i))
    if show: plt.show()
    
def test0(filepath):
    print('-----------------------------------------Test0---')
    lattice = factory(filepath)
    # dictprnt(PARAMS, 'PARAMS')
    sollTrack = track_soll(lattice)
    DEBUG_TEST0('sollTrack:\n'+sollTrack.asTable())
    DEBUG_TEST0('sollTrack:\n(first): {}\n (last): {}'.format(sollTrack.first_str(),sollTrack.last_str()))

# def test1(filepath):
def test1(options):
    print('-----------------------------------------Test1---')
    tracker(options)
    
if __name__ == '__main__':
    DEBUG_TRACK       = DEBUG_OFF
    DEBUG_SOLL_TRACK  = DEBUG_OFF
    DEBUG_TEST0       = DEBUG_ON
    
    options = dict( filepath = 'yml/work.yml',
                    particlesPerBunch = 3000,
                    show    = True,
                    save    = False,
                    skip    = 4
                    )
    template = Template('$tx1 $tx2 $tx3 $tx4')
    # filepath = 'yml/work.yml'
    # test0(filepath)
    test1(options)
