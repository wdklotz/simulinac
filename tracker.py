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
from setutil import WConverter
from bunch import BunchFactory, Gauss1D, Track, Tpoint
from trackPlot import poincarePlot

# DEBUGGING
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass

#todo: show invalid tracks
def scatterPlot(bunch, ordinate, abszisse, text, minmax=(1.,1.)):
    """ prepare the plot of a Poincaré section """
    
    txt = ('{} initial'.format(text),'{} final'.format(text))
    initial = 0
    final   = -1
    for loc in (initial,final): # figures
        x=[]; y=[]
        fig = plt.figure(loc)
        box = '{} {} particles'.format(txt[loc],bunch.nbofparticles())
        for nd, particle in iter(bunch): # particles
            track = particle['track']
            nd, tpoint = track.getpoints()[loc]
            point = tpoint()
            x.append(point[ordinate])
            y.append(point[abszisse])
        xmax = max(x)
        xmin = min(x)
        ymax = max(y)
        ymin = min(y)
        minmax = (xmax-xmin,ymax-ymin)
        # plt.scatter(x,y,s=1)
        poincarePlot((x,y), box, max = minmax, projections = (1,1))
    return
    
def projection(bunch, ordinate = K.z, abszisse = K.zp, show = True, save = False):
    """ 
    2D phase space projections of Poincaré sections 
    """
    symbol    = ("x","x'","y","y'","z","z'")
#todo: correct calc of minmax
    # sigmas    = (bunch['sigx'],bunch['sigxp'],bunch['sigy'],bunch['sigyp'],bunch['sigz'],bunch['sigzp'])
    text      = '{}-{}'.format(symbol[ordinate],symbol[abszisse])
    # minmax    = (10.e3*sigmas[ordinate],10.e3*sigmas[abszisse])
    
    fig = scatterPlot(bunch, ordinate=ordinate, abszisse=abszisse, text=text)
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
    nb,last_point = track.getpoints()[-1]
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
        pass
    for node in iter(lattice):              # nodes
        ncnt +=1
        for nd, particle in iter(bunch):    # particles
            if not particle.lost: 
                failure = track_node(node,particle)
                if failure:
                    particle.lost = True
                    lost += 1
                    continue
                # showing track-loop progress
                if (ncnt+1)%lmod == 0:
                    tx4    = ' {}% done {}/{} lost/initial'.format(int(ncnt/lnode*100), lost, nbpart)
                    # tx = ('(track design)', '(track bunch)', zeuge[ncnt%4], tx4)
                    tx = ('(track design)', '(track bunch)', '', tx4)
                    progress(tx)
            else:
                continue
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
    npart             = options['particles_per_bunch']
    show              = options['show']
    save              = options['save']
    skip              = options['skip']
    tkin              = PARAMS['sollteilchen'].tkin
    conv              = WConverter(tkin)

    #make lattice
    t0 = time.clock()
    lattice = factory(filepath)
    # calculate twiss paramters at entrance
    waccept(lattice.first_gap)
    t1 = time.clock()

    # bunch-configuration from PARAMS
    # {x(x)xp}
    twx = PARAMS['twiss_x_i']
    betax_i,alfax_i,gammax_i,emitx_i = twx()
    sigma_x   = twx.sigmaH()
    sigma_xp  = twx.sigmaV()
    # {y(x)yp}
    twy = PARAMS['twiss_y_i']
    betay_i,alfay_i,gammay_i,emity_i = twy()
    sigma_y   = twy.sigmaH()
    sigma_yp  = twy.sigmaV()
    # {Dphi(x)w}
    tww = PARAMS['twiss_w_i']
    betaw,alfaw,gammaw,emitw = tww()
    sigma_Dphi  = tww.sigmaH()
    sigma_w     = tww.sigmaV()
    # {z(x)Dp2p}
    twz = PARAMS['twiss_z_i']
    betaz,alfaz,gammaz,emitz = twz()
    sigma_z    = twz.sigmaH()
    sigma_Dp2p = twz.sigmaV()
#todo: working
    # gather for print
    parameter_log = {}
    parameter_log['tkin [MeV]'] = tkin
    parameter_log["sigma(x,x')_i"] = (sigma_x,sigma_xp)
    parameter_log["sigma(y,y')_i"] = (sigma_y,sigma_yp)
    parameter_log["sigma(Dphi,w)_i"] = (sigma_Dphi,sigma_w)
    parameter_log['betax_i'] = betax_i
    parameter_log['betay_i'] = betay_i
    parameter_log['betaw_i'] = betaw
    parameter_log['emitx_i'] = emitx_i
    parameter_log['emity_i'] = emity_i
    parameter_log['emitw_i'] = emitw
    dictprnt(parameter_log,'Tracker Options'); print()

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

# def track_bunch(options):
#     print('-----------------------------------------track_bunch---')
#     tracker(options)
    
if __name__ == '__main__':
    DEBUG_TRACK       = DEBUG_OFF
    DEBUG_SOLL_TRACK  = DEBUG_OFF
    DEBUG_TEST0       = DEBUG_ON
    
    template_file = 'yml/worktmpl.yml'          # template for m4
    input_file    = 'yml/trackIN.yml'           # input for tracker.py
    macros_file   = 'yml/macros.sh'              # macro definitions
    command = "{} {} > {}".format(macros_file,template_file, input_file)
    print('m4->script: ',macros_file,' template: ',template_file,' input: ',input_file)
    os.system(command)

    # test0(input_file)

    options = dict( input_file = input_file,
                    particles_per_bunch = 1000,
                    show    = True,
                    save    = False,
                    skip    = 1
                    )
    tracker(options)
