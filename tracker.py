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
# import matplotlib.mlab as mlab
# import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import time
from trackPlot import poincarePlot
from string import Template

from lattice_generator import Factory
import elements as ELM
from setutil import DEBUG, PARAMS, dictprnt, sigmas, K
from bunch import Track, Bunch, Gauss1D

# DEBUGGING
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_TRACK       = DEBUG_OFF
DEBUG_SOLL_TRACK  = DEBUG_OFF
DEBUG_TEST0       = DEBUG_ON

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
    boxtext = '{} {} particles'.format(txt, bunch.nbtracks)
    poincarePlot(x, y, boxtext, minmax, projections = (1,1))
    return fig

progress = Template('$tx1 $tx2 $tx3')

def track(lattice,bunch):
    """
    Tracks a bunch of particles through the lattice using maps
    - lattice is a list of elements
    - bunch   is a list of particle tracks
    """
    from time import sleep

    if DEBUG_TRACK == DEBUG_ON: dictprnt(bunch._params,'bunch',filter='tracklist')
    
    invalid_tracks = []
    valid_tracks   = []
    losses         = 0

    for (tcount, ptrack) in enumerate(bunch.tracks):       # loop tracks
        invalid = False
        ti = ptrack.first()

        for element in lattice.seq:                        # loop elements
            try:
                tf = element.map(ti)                       # map!
            except ValueError as ex:
                invalid_tracks.append(ptrack)
                losses = len(invalid_tracks)
                invalid = True
                break
            ti = tf
            # take a poincare section at MRK elements only
            if isinstance(element,ELM.MRK):
                ptrack.append(tf)

        if not invalid : valid_tracks.append(ptrack)
        # showing some track-loop cycles progress
        if (tcount+1)%25 == 0:
            prog = progress.substitute(tx1='(soll-track)', tx2='(tracks)', tx3='{}/{}/{} done/lost/initial'.format(tcount+1,losses,bunch.nbtracks))
            print('\r{}'.format(prog), end='')

    # keep valid tracks in the bunch
    bunch.tracks = valid_tracks
    return bunch

def track_soll(lattice):
    """
    Tracks the reference particle through the lattice and redefines the lattice element parameters to
    adapted to the energy of the accellerated reference particle.
    """
    soll_track = Track(start=np.array([ 0., 0., 0., 0., 0., 0., PARAMS['sollteilchen'].tkin, 1., 0., 1.]))
    for element in lattice.seq:
        DEBUG_SOLL_TRACK(element,' pos {:.4f} label "{}"'.format(element.position[1],element.label))
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

def tracker(filepath, particlesPerBunch, show, save, skip):
    """ prepare and launch tracking """
    #make lattice with time
    t0 = time.clock()
    lattice = Factory(filepath)
    t1 = time.clock()

    # bunch-configuration
    alfax_i   = PARAMS['alfax_i'] 
    betax_i   = PARAMS['betax_i']
    emitx_i   = PARAMS['emitx_i']
    alfay_i   = PARAMS['alfay_i']
    betay_i   = PARAMS['betay_i']
    emity_i   = PARAMS['emity_i']
    sigmaz_i  = PARAMS['sigmaz_i']
    sigmazp_i = PARAMS['dp2p_i']

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
    prog = progress.substitute(tx1='(soll-track)', tx2='', tx3='')
    print('\r{}'.format(prog), end='')
    t2 = time.clock()
    track_soll(lattice)  # track soll
    t3 = time.clock()
    prog = progress.substitute(tx1='(soll-track)', tx2='(tracks)', tx3='')
    print('\r{}'.format(prog), end='')
    track(lattice,bunch) # track bunch
    t4 = time.clock()
    # make 2D projections
    # track_plane(bunch, K.x, K.xp, show, save)
    # track_plane(bunch, K.x, K.y, show, save)
    track_plane(bunch, K.z, K.zp, show, save, skip)
    t5 = time.clock()

    print()
    print('total time     >> {:6.3f} [sec]'.format(t5-t0))
    print('parse lattice  >> {:6.3f} [sec] {:4.1f} [%]'.format((t1-t0),(t1-t0)/(t5-t0)*1.e2))
    print('generate bunch >> {:6.3f} [sec] {:4.1f} [%]'.format((t2-t1),(t2-t1)/(t5-t0)*1.e2))
    print('track design   >> {:6.3f} [sec] {:4.1f} [%]'.format((t3-t2),(t3-t2)/(t5-t0)*1.e2))
    print('track bunch    >> {:6.3f} [sec] {:4.1f} [%]'.format((t4-t3),(t4-t3)/(t5-t0)*1.e2))
    print('fill plots     >> {:6.3f} [sec] {:4.1f} [%]'.format((t5-t4),(t5-t4)/(t5-t0)*1.e2))

def track_plane(bunch, ordinate, abzisse, show, save, skip):
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
    lattice = Factory(filepath)
    # dictprnt(PARAMS, 'PARAMS')
    sollTrack = track_soll(lattice)
    DEBUG_TEST0('sollTrack:\n'+sollTrack.asTable())
    DEBUG_TEST0('sollTrack:\n(first): {}\n (last): {}'.format(sollTrack.first_str(),sollTrack.last_str()))

def test1(filepath):
    print('-----------------------------------------Test1---')
    print('tracker() with lattice-file {}'.format(filepath))
    tracker(filepath, particlesPerBunch = 3000, show=False, save=True, skip=1)
    
if __name__ == '__main__':
    filepath = 'yml/work.yml'    ## the default input file (YAML syntax)
    # test0(filepath)
    test1(filepath)
