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
from setutil import DEBUG, PARAMS, dictprnt, sigmas, K, PARAMS, waccept
from bunch import Track, Bunch, Gauss1D, Tpoint
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
    nbmkrs = bunch.nbmkrs        # poincare markers per track
    psec = poincare_section
    if psec == 0:
        txt = '{} initial'.format(text)
    elif psec == (nbmkrs-1):
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
    
# def track_track(track,lattice):
#     """ successive mappings of a track through the lattice """
#     invalid = False
# 
#     ti = track.first()
#     for element in lattice.seq:                  # loop elements
#         try:
#             tf = element.map(ti)                 # map!
#         except ValueError as ex:
#             invalid = True
#             break
#         ti = tf
#         # take a poincare section at MRK elements with poincare action only
#         if isinstance(element,ELM.MRK) and not invalid:
#             if 'poincare' in element.actions:
#                 track.append(tf)
#     return (track,invalid)

def track_node(node,particle):
    """
    Tracks a particle through a node
    """
    lost = False
    track = particle['track']
    n,last_point = track.getpoints()[-1]
    try:
        new_point = node.map(last_point)
    except ValueError as ex:
        lost = True
    track.addpoint(new_point)
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
    # if DEBUG_TRACK == DEBUG_ON: dictprnt(bunch._params,'bunch',filter='tracklist')
    
    losses         = 0
    zeuge          = ('*\u007C*','**\u007C','*\u007C*','\u007C**')  # *|*
    tx4            = ' {}/{}/{} done/lost/initial'.format(0,0,bunch.nbtracks)

    lost_particles = []
    live_particles = []
    
    for node in iter(lattice):
        pcount = 0
        for pnumber, particle in iter(bunch):    # loop particles in bunch
            pcount += 1
            lost = track_node(node,particle)
            if lost:
                lost_particles.append(pnumber)
            else:
                live_particles.append(pnumber)
            # showing track-loop progress
            if (pcount+1)%25 == 0:
                losses = len(lost_particles)
                tx4    = ' {}/{}/{} done/lost/initial'.format(tcount+1, losses, bunch.nbtracks)
            progress(('(track design)', '(track bunch)', zeuge[tcount%4], tx4))
    # log valid/invalid tracks in the bunch
    bunch.setlost(lost_particles)
    bunch.setlive(live_particles)
    print('\ndone: live particles {}, lost paricles {}'.format(len(live_particles),len(lost_particles)))
    return bunch

def track_soll(lattice):
    """
    Tracks the reference particle through the lattice and redefines the lattice element parameters to
    adapted to the energy of the accelerated reference particle.
    """
    soll_track = Track()
    soll_track.addpoint(NP.array([ 0., 0., 0., 0., 0., 0., PARAMS['sollteilchen'].tkin, 1., 0., 1.]))
    for node in iter(lattice):
        # DEBUG_SOLL_TRACK(element,' pos {:.4f} label "{}"'.format(element.position[1],element.label))
        n,pi = soll_track.getpoints()[-1]     # point: at entrance
        # DEBUG_SOLL_TRACK('track_soll(i) ', soll_track.point_str(ti))
        # DEBUG_SOLL_TRACK('track_soll: complete track\n{}'.format(soll_track.points_str()))
        """ energy adjustment """
        element.adjust_energy(pi[K.T])
        # DEBUG_SOLL_TRACK('track_soll: complete track\n{}'.format(soll_track.points_str()))
        """ mapping with soll map """
        pf = element.soll_map(pi)
        # DEBUG_SOLL_TRACK('track_soll(f) ',soll_track.point_str(tf))
        soll_track.addpoint(pf)
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

    #make lattice
    t0 = time.clock()
    lattice = factory(filepath)
    t1 = time.clock()

    # bunch-configuration from PARAMS
    alfax_i   = PARAMS['alfax_i'] 
    betax_i   = PARAMS['betax_i']
    emitx_i   = PARAMS['emitx_i']
    alfay_i   = PARAMS['alfay_i']
    betay_i   = PARAMS['betay_i']
    emity_i   = PARAMS['emity_i']
    # calculate longitudinal paramters at entrance
    waccept(lattice.first_gap)
    alfaz_i   = 0.
    betaz_i   = PARAMS['betaz_i']
    emitz_i   = PARAMS['emitz_i']

    sigmas_x  = sigmas(alfax_i, betax_i, emitx_i)
    sigmas_y  = sigmas(alfay_i, betay_i, emity_i)
    sigmas_z  = sigmas(      0.,betaz_i, emitz_i)

    options['tkin [MeV'] = PARAMS['sollteilchen'].tkin
    options["sigma(x,x')_i"] = sigmas_x
    options["sigma(y,y')_i"] = sigmas_y
    options["sigma(z,z')_i"] = sigmas_z
    dictprnt(options,'Tracker Options'); print()

#todo: need new bunch factory
    bunch = Bunch()
    bunch.nbtracks    = particlesPerBunch
    bunch.coord_mask  = (1,1,1,1,1,1)
    bunch.disttype    = Gauss1D
    bunch['sigma-x']  = sigmas_x[0]
    bunch['sigma-xp'] = sigmas_x[1]
    bunch['sigma-y']  = sigmas_y[0]
    bunch['sigma-yp'] = sigmas_y[1]
    bunch['sigma-z']  = sigmas_z[0]
    bunch['sigma-zp'] = sigmas_z[1]
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
    nbmkrs  = [i for i in range(bunch.nbmkrs)]
    # print('nbmkrs',nbmkrs)
    for i in nbmkrs:
        if i != 0 and i != nbmkrs[-1] and i%skip != 0: continue   # skip nbmkrs, keep first and last
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
def track_bunch(options):
    print('-----------------------------------------Test1---')
    tracker(options)
    
if __name__ == '__main__':
    DEBUG_TRACK       = DEBUG_OFF
    DEBUG_SOLL_TRACK  = DEBUG_OFF
    DEBUG_TEST0       = DEBUG_ON
    
    options = dict( filepath = 'yml/work.yml',
                    particlesPerBunch = 1000,
                    show    = True,
                    save    = False,
                    skip    = 1
                    )
    template = Template('$tx1 $tx2 $tx3 $tx4')
    # filepath = 'yml/work.yml'
    # test0(filepath)
    track_bunch(options)
