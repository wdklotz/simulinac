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
import os
import numpy as NP
import matplotlib.pyplot as plt
import time
from string import Template
from math import sqrt

from lattice_generator import factory
import elements as ELM
from setutil import DEBUG, PARAMS, dictprnt, sigmas, K, PARAMS, waccept
from bunch import BunchFactory, Gauss1D, Track, Tpoint
from trackPlot import poincarePlot

# DEBUGGING
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass

def scatterPlot(bunch, ordinate, abszisse, text, minmax):
    """ prepare the plot of a Poincaré section """
    txt = '{} final'.format(text)
    x=[]; y=[]
    fig = plt.figure(-1)
    for nparticle,particle in iter(bunch):
        track = particle['track']
        ntp, tpoint = track[-1]
        point = tpoint()
        x.append(point[ordinate]*1.e3)
        y.append(point[abszisse]*1.e3)

#todo: invalid tracks
    # xi=[]; yi=[]
    # for t in bunch.invalid_tracks:             # loop invalid tracks
    #     if psec < t.nbpoints:
    #         xi.append(t.point_at(psec)[ordinate]*1.e3)
    #         yi.append(t.point_at(psec)[abszisse]*1.e3)
    # bad = (xi,yi)

    box = '{} {} particles'.format(txt, bunch.nbofparticles())
    poincarePlot((x,y), box, max = minmax, projections = (1,1))
    return fig
    
def projection(bunch, ordinate = K.z, abszisse = K.zp, show = True, save = False):
    """ 
    2D phase space projections of Poincaré sections 
    """
    symbol    = ("x","x'","y","y'","z","z'")
    sigmas    = (bunch['sigx'],bunch['sigxp'],bunch['sigy'],bunch['sigyp'],bunch['sigz'],bunch['sigzp'])
    text      = '{}-{}'.format(symbol[ordinate],symbol[abszisse])
    minmax    = (10.e3*sigmas[ordinate],10.e3*sigmas[abszisse])
    
    fig = scatterPlot(bunch, ordinate=ordinate, abszisse=abszisse, text=text, minmax=minmax)
    if save: plt.savefig('figures/poincare_section_{}_{}.png'.format(text,i))
    if show: plt.show()

def progress(tx):
    template = Template('$tx1 $tx2 $tx3 $tx4')
    res = template.substitute(tx1=tx[0] , tx2=tx[1] , tx3=tx[2] , tx4=tx[3] )
    print('{}\r'.format(res),end='')

#todo: needs aperture check here
def track_node(node,particle):
    """
    Tracks a particle through a node
    """
    lost = False
    track = particle['track']
    n,last_point = track.getpoints()[-1]
    try:
        new_point = node.map(last_point())
    except ValueError as ex:
        lost = True
    track.addpoint(Tpoint(point=new_point))
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

    ncnt   = 0
    lnode  = len(lattice.seq)
    lmod   = int(lnode*0.05)
    lost   = 0
    nbpart = bunch.nbofparticles()
    for node in iter(lattice):
        ncnt += 1
        for pnumber, particle in iter(bunch):    # loop particles in bunch
            if particle.lost: continue
            if track_node(node,particle):
                particle.lost = True
                lost += 1
            # showing track-loop progress
            if (ncnt+1)%lmod == 0:
                tx4    = ' {}% done {}/{} lost/initial'.format(int(ncnt/lnode*100), lost, nbpart)
                # tx = ('(track design)', '(track bunch)', zeuge[ncnt%4], tx4)
                tx = ('(track design)', '(track bunch)', '', tx4)
                progress(tx)
    live = nbpart - lost
    print('\ndone: live particles {}, lost particles {}'.format(live,lost))
    return bunch

def track_soll(lattice):
    """
    Tracks the reference particle through the lattice and redefines the lattice 
    element parameters according to the energy of the accelerated reference particle.
    """
    soll_track = Track()
    tpoint = Tpoint(NP.array([ 0., 0., 0., 0., 0., 0., PARAMS['sollteilchen'].tkin, 1., 0., 1.]))
    soll_track.addpoint(tpoint)
    for node in iter(lattice):
        n,pi = soll_track.getpoints()[-1]   # n,Tpoint at entrance
        """ energy adjustment """
        node.adjust_energy(pi()[K.T])
        """ mapping with soll map """
        pf = node.soll_map(pi())
        tpoint = Tpoint(pf)                  # n+1,Tpoint at exit
        soll_track.addpoint(tpoint)
    return soll_track

def tracker(options):
    """ prepare and launch tracking """
    filepath          = options['input_file']
    particles_per_bunch = options['particles_per_bunch']
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
    options['tkin [MeV]'] = PARAMS['sollteilchen'].tkin
    options["sigma(x,x')_i"] = sigmas_x
    options["sigma(y,y')_i"] = sigmas_y
    options["sigma(z,z')_i"] = sigmas_z
    dictprnt(options,'Tracker Options'); print()
    # bunch factory
    bunchfactory = BunchFactory(Gauss1D)
    bunchfactory['sigx']        = sigmas_x[0]
    bunchfactory['sigxp']       = sigmas_x[1]
    bunchfactory['sigy']        = sigmas_y[0]
    bunchfactory['sigyp']       = sigmas_y[1]
    bunchfactory['sigz']        = sigmas_z[0]
    bunchfactory['sigzp']       = sigmas_z[1]
    bunchfactory['coord_mask']  = (1,1,1,1,1,1)
    bunchfactory['nbparticles'] = particles_per_bunch
    bunchfactory['tkin']        = PARAMS['injection_energy']
    bunch = bunchfactory()
    # launch tracking and show final with time
    progress(('(track design)', '', '', ''))
    t2 = time.clock()
    track_soll(lattice)  # <----- track soll
    t3 = time.clock()
    progress(('(track design)', '(track bunch)', '', ''))
    track(lattice,bunch) # <----- track bunch
    t4 = time.clock()
    # make 2D projections
    projection(bunch, ordinate = K.z, abszisse = K.zp, show = True, save = False)
    t5 = time.clock()
    # finish up
    print()
    print('total time     >> {:6.3f} [sec]'.format(t5-t0))
    print('parse lattice  >> {:6.3f} [sec] {:4.1f} [%]'.format((t1-t0),(t1-t0)/(t5-t0)*1.e2))
    print('generate bunch >> {:6.3f} [sec] {:4.1f} [%]'.format((t2-t1),(t2-t1)/(t5-t0)*1.e2))
    print('track design   >> {:6.3f} [sec] {:4.1f} [%]'.format((t3-t2),(t3-t2)/(t5-t0)*1.e2))
    print('track bunch    >> {:6.3f} [sec] {:4.1f} [%]'.format((t4-t3),(t4-t3)/(t5-t0)*1.e2))
    print('fill plots     >> {:6.3f} [sec] {:4.1f} [%]'.format((t5-t4),(t5-t4)/(t5-t0)*1.e2))

def test0(filepath):
    print('-----------------------------------------Test0---')
    lattice = factory(filepath)
    sollTrack = track_soll(lattice)
    table = sollTrack.as_table()
    DEBUG_TEST0('sollTrack:\n'+table)
    d,first = sollTrack[0]
    d,last  = sollTrack[-1]
    DEBUG_TEST0('sollTrack:\n(first): {}\n (last): {}'.format(first.as_str(),last.as_str()))

# def test1(filepath):
def track_bunch(options):
    print('-----------------------------------------track_bunch---')
    tracker(options)
    
if __name__ == '__main__':
    DEBUG_TRACK       = DEBUG_OFF
    DEBUG_SOLL_TRACK  = DEBUG_OFF
    DEBUG_TEST0       = DEBUG_ON
    
    template_file = 'yml/worktmpl.yml'
    input_file    = 'yml/trackIN.yml'
    preproc_file  = 'yml/ppdef.sh'
    command = "{} {} > {}".format(preproc_file,template_file, input_file)
    print('m4->script: ',preproc_file,' template: ',template_file,' input: ',input_file)
    os.system(command)

    # test0(input_file)

    options = dict( input_file = input_file,
                    particles_per_bunch = 3000,
                    show    = True,
                    save    = False,
                    skip    = 1
                    )
    track_bunch(options)
