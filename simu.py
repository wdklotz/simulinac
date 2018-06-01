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

    You should have received a copy of the GNU General Public Licensedir
    along with SIMULINAC.  If not, see <http://www.gnu.org/licenses/>.
"""
#todo: update simu_manual.odt
#todo: update README.md
#todo: revise flow control with FLAGS (partly done)
import sys
from math import sqrt
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from setutil import PARAMS,FLAGS,SUMMARY,dictprnt,DEBUG
from setutil import collect_data_for_summary, waccept, elli
from lattice_generator import parse_yaml_and_fabric
from tracker import track_soll

import bucket_size

# DEBUG
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_MODULE  = DEBUG_OFF

def display(functions):
    figures = []
    plots   = []
    if FLAGS['csTrak'] and FLAGS['dWf'] == 0:
        plots.append(display0) # CS tracks {x,y}
    elif FLAGS['csTrak'] and FLAGS['dWf'] == 1:
        plots.append(display1) # CS tracks {x,y,z}
        if FLAGS['bucket']:
            plots.append(bucket) # separatrix
    if FLAGS['ellipse']:
        plots.append(display2) # ellipses

    # all plots and their figures
    if len(plots) != 0:
        print('PREPARE DISPLAY')
        [plot(figures,functions) for plot in plots]
        [plt.show(fig) for fig in figures]

def bucket(figures,dummy):
    bucket_size.bucket(figures)
    
def display0(figures,functions):
    """
    CS-Tracks w/o longitudinal motion
    """
    #----------*----------*   # unpack
    sigm_fun = functions[0]
    cos_like = functions[1]
    sin_like = functions[2]
    lat_plot = functions[3]
    #-------------------- Bahnkoordinate (z)
    z    = [x[0] for x in sigm_fun]    # Ordinate
    bx   = [x[1] for x in sigm_fun]    # envelope (sigma-x)
    by   = [x[2] for x in sigm_fun]    # envelope (sigma-y)
    zero = [0.   for x in sigm_fun]    # zero line
    #-------------------- trajectories (tz)
    tz  = [x[0] for x in cos_like]   # Ordinate
    cx  = [x[1] for x in cos_like]   # cos-like-x
    cy  = [x[3] for x in cos_like]   # cos-like-y
    cz  = [x[5] for x in cos_like]   # cos-like-z
    cdw = [x[6] for x in cos_like]   # cos-like-dw/w
    sx  = [x[1] for x in sin_like]   # sin-like-x
    sy  = [x[3] for x in sin_like]   # sin-like-x
    sz  = [x[5] for x in sin_like]   # sin-like-z
    sdw = [x[6] for x in sin_like]   # sin-like-dw/w
    #-------------------- lattice viseo
    stop_viseo = 2000                  # stop viseo plot after so many points
    vis_ordinate = [x[0] for x in lat_plot]
    vis_abszisse = [x[1] for x in lat_plot]
    vzero        = [0.   for x in lat_plot]      # zero line
    #-------------------- figure frame
    width=14; height=7.6
    # fighdr = 'lattice version = {}, input file = {}'.format(PARAMS['lattice_version'],PARAMS['input_file'])
    fig = plt.figure(num=0,figsize=(width,height),facecolor='#eaecef',tight_layout=True)
    #-------------------- transverse X
    splot=plt.subplot(211)
    splot.set_title('transverse x')
    plt.plot(z,bx ,label=r'$\sigma$ [m]',color='green')
    plt.plot(tz,cx,label='Cx[m]',color='blue',linestyle='-')
    plt.plot(tz,sx,label='Sx[m]',color='red' ,linestyle='-')
    vscale=plt.axis()[3]*0.1
    viseox = [x*vscale for x in vis_abszisse]
    for i in range(stop_viseo,len(vis_ordinate)): viseox[i] = 0.   # stop lattice plotting
    plt.plot(vis_ordinate,viseox,label='',color='black')
    plt.plot(vis_ordinate,vzero,color='black')
    plt.legend(loc='lower right',fontsize='x-small')
    #-------------------- transverse Y
    splot=plt.subplot(212)
    splot.set_title('transverse y')
    plt.plot(z,by ,label=r'$\sigma$ [m]',color='green')
    plt.plot(tz,cy,label='Cy[m]',color='blue',linestyle='-')
    plt.plot(tz,sy,label='Sy[m]',color='red' ,linestyle='-')
    vscale=plt.axis()[3]*0.1
    viseoy = [x*vscale for x in vis_abszisse]
    for i in range(stop_viseo,len(vis_ordinate)): viseoy[i] = 0.   # stop lattice plotting
    plt.plot(vis_ordinate,viseoy,label='',color='black')
    plt.plot(vis_ordinate,vzero,color='black')
    plt.legend(loc='lower right',fontsize='x-small')

    figures.append(fig)

def display1(figures,functions):
    """
    CS-Tracks with longitudinal motion
    """
    #-------------------- unpack
    sigm_fun = functions[0]
    cos_like = functions[1]
    sin_like = functions[2]
    lat_plot = functions[3]
    #-------------------- twiss functions
    z    = [x[0] for x in sigm_fun]  # Ordinate
    bx   = [x[1] for x in sigm_fun]  # envelope (sigma-x)
    by   = [x[2] for x in sigm_fun]  # envelope (sigma-y)
    zero = [0.   for x in sigm_fun]  # zero line
    #-------------------- trajectories
    tz  = [x[0] for x in cos_like]   # Ordinate
    cx  = [x[1] for x in cos_like]   # cos-like-x
    cy  = [x[3] for x in cos_like]   # cos-like-y
    cz  = [x[5] for x in cos_like]   # cos-like-z
    cdw = [x[6] for x in cos_like]   # cos-like-dw/w
    sx  = [x[1] for x in sin_like]   # sin-like-x
    sy  = [x[3] for x in sin_like]   # sin-like-x
    sz  = [x[5] for x in sin_like]   # sin-like-z
    sdw = [x[6] for x in sin_like]   # sin-like-dw/w
    #-------------------- lattice viseo
    stop_viseo = 2000                  # stop viseo plot after so many points
    vis_ordinate = [x[0] for x in lat_plot]
    vis_abszisse = [x[1] for x in lat_plot]
    vzero        = [0.   for x in lat_plot]      # zero line
    #-------------------- figure frame
    width=14; height=7.6
    # fighdr = 'lattice version = {}, input file = {}'.format(PARAMS['lattice_version'],PARAMS['input_file'])
    fig = plt.figure(num=1,figsize=(width,height),facecolor='#eaecef',tight_layout=True)
    #-------------------- transverse X
    splot=plt.subplot(311)
    splot.set_title('transverse x')
    plt.plot(z,bx ,label=r'$\sigma$ [m]',color='green')
    plt.plot(tz,cx,label='C [m]',color='blue',linestyle=':')
    plt.plot(tz,sx,label='S [m]',color='red' ,linestyle=':')
    vscale=plt.axis()[3]*0.4
    viseox = [x*vscale for x in vis_abszisse]
    for i in range(stop_viseo,len(vis_ordinate)): viseox[i] = 0.   # stop lattice plotting
    plt.plot(vis_ordinate,viseox,label='',color='black')
    plt.plot(vis_ordinate,vzero,color='black')
    plt.legend(loc='lower right',fontsize='x-small')
    #-------------------- transverse Y
    splot=plt.subplot(312)
    splot.set_title('transverse y')
    plt.plot(z,by ,label=r'$\sigma$ [m]',color='green')
    plt.plot(tz,cy,label='C [m]',color='blue',linestyle=':')
    plt.plot(tz,sy,label='S [m]',color='red' ,linestyle=':')
    vscale=plt.axis()[3]*0.4
    viseoy = [x*vscale for x in vis_abszisse]
    for i in range(stop_viseo,len(vis_ordinate)): viseoy[i] = 0.   # stop lattice plotting
    plt.plot(vis_ordinate,viseoy,label='',color='black')
    plt.plot(vis_ordinate,vzero,color='black')
    plt.legend(loc='lower right',fontsize='x-small')
    #-------------------- longitudinal dPhi, dW/W
    ax_l=plt.subplot(313)
    ax_l.set_title('longitudinal')
    ax_l.set_ylabel(r"$\Delta\phi$ [deg]")
    ax_l.tick_params(axis='y', colors='green')
    ax_l.yaxis.label.set_color('green')
    ax_l.plot(tz,cz,color='green',linestyle=':')
    # ax_l.plot(tz,sz,color='green')
    vscale=ax_l.axis()[3]*0.7
    viseoz = [x*vscale for x in vis_abszisse]
    for i in range(stop_viseo,len(vis_ordinate)): viseoz[i] = 0.   # stop lattice plotting
    ax_l.plot(vis_ordinate[0:stop_viseo],viseoz[0:stop_viseo],label='',color='black')
    ax_l.plot(vis_ordinate,vzero,color='green',linestyle='--')

    ax_r = ax_l.twinx()
    ax_r.set_ylabel(r'$\Delta$w/w [%]')
    ax_r.tick_params(axis='y', colors='red')
    ax_r.yaxis.label.set_color('red')
    # ax_r.plot(tz,cdw,color='red',linestyle=':')
    ax_r.plot(tz,sdw,color='red')
    ax_r.plot(vis_ordinate,vzero,color='red', linestyle='--')

    figures.append(fig)

def display2(figures,dummy):
    """ display x- and y-phase-space ellipses """
    xy = (0,0)
    ellix = elli(xy,PARAMS['alfax_i'],PARAMS['betax_i'],PARAMS['emitx_i'])
    elliy = elli(xy,PARAMS['alfay_i'],PARAMS['betay_i'],PARAMS['emity_i'])

    ells = [Ellipse(*ellix,color='green',fill=False),Ellipse(*elliy,color='red',fill=False)]

    fig, ax = plt.subplots()
    fig.suptitle('transverse phase-space {[m],[rad]}')
    fig.legend(ells,("{x,x'}","{y,y'}"),loc=1)

    for e in ells:
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)

    x1 = sqrt(PARAMS['emitx_i']*PARAMS['betax_i'])
    x2 = sqrt(PARAMS['emity_i']*PARAMS['betay_i'])
    xmax = max(x1,x2)
    gammax = (1.+PARAMS['alfax_i']**2)/PARAMS['betax_i']
    gammay = (1.+PARAMS['alfay_i']**2)/PARAMS['betay_i']
    y1 = sqrt(PARAMS['emitx_i']*gammax)
    y2 = sqrt(PARAMS['emity_i']*gammay)
    ymax = max(y1,y2)
    scale = 0.6
    plt.xlim(-xmax*scale, xmax*scale)
    plt.ylim(-ymax*scale, ymax*scale)    

    figures.append(fig)
        
#                            |----------------------- |
# -------------------------  | everything starts here |
#                            |----------------------- |
def simulation(filepath):
    # parse input file and create a lattice
    lattice = parse_yaml_and_fabric(filepath)

    # calculate longitudinal paramters at entrance
    waccept(lattice.first_gap)

    # configure elements with increasing energy
    soll_track = track_soll(lattice)

    # count elements and make other statistics
    lattice.stats(soll_track)

    # ganzer Beschleuniger, Anfangswerte, Summary, etc...
    lattice.cell(closed = FLAGS['periodic'])

    # make summary
    collect_data_for_summary(lattice)

    # grafics & results
    lat_plot = lattice.lattice_plot_function()
    kv_only  = FLAGS['KVout']
    steps    = 23
    if kv_only: 
        dictprnt(PARAMS,text='PARAMS',njust=1)
    else:
        (c_like,s_like) = lattice.cs_traj(steps = steps) # calc sin- and cos-like trajectories
        sigma_fun       = lattice.sigmas(steps = steps)
        functions       = (sigma_fun,c_like,s_like,lat_plot)
        # show summary
        dictprnt(SUMMARY,text='summary')
        # call plot dispatcher
        display(functions)

if __name__ == '__main__':
    filepath = 'yml/ref_run.yml'# the default input file (YAML syntax)
    filepath = 'yml/work.yml'

    if len(sys.argv) == 2:
        filepath = sys.argv[1]
    simulation(filepath)
    
