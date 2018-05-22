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
import sys
from math import sqrt
from matplotlib.pyplot import plot,show,legend,figure,subplot,axis

from setutil import PARAMS,FLAGS,SUMMARY,dictprnt,DEBUG
from setutil import collect_data_for_summary
from lattice_generator import parse_yaml_and_fabric
from tracker import track_soll

# DEBUG
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_MODULE = DEBUG_OFF
DEBUG_BUCKET = DEBUG_OFF

def display(functions):
    print('PREPARE DISPLAY')
    if FLAGS['dWf'] == 0:
        display0(functions)
    else:
        display1(functions)
        if DEBUG_BUCKET == DEBUG_ON:
            from bucket_size import bucket
            bucket()     # separatrix

def display0(functions):
    """
    Plotting w/o longitudinal motion
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
    fighdr = 'lattice version = {}, input file = {}'.format(PARAMS['lattice_version'],PARAMS['input_file'])
    figure(fighdr,figsize=(width,height),facecolor='#eaecef',tight_layout=True)
    #-------------------- transverse X
    splot=subplot(211)
    splot.set_title('transverse x')
    plot(z,bx ,label=r'$\sigma$ [m]',color='green')
    plot(tz,cx,label='Cx[m]',color='blue',linestyle='-')
    plot(tz,sx,label='Sx[m]',color='red' ,linestyle='-')
    vscale=axis()[3]*0.1
    viseox = [x*vscale for x in vis_abszisse]
    for i in range(stop_viseo,len(vis_ordinate)): viseox[i] = 0.   # stop lattice plotting
    plot(vis_ordinate,viseox,label='',color='black')
    plot(vis_ordinate,vzero,color='black')
    legend(loc='lower right',fontsize='x-small')
    #-------------------- transverse Y
    splot=subplot(212)
    splot.set_title('transverse y')
    plot(z,by ,label=r'$\sigma$ [m]',color='green')
    plot(tz,cy,label='Cy[m]',color='blue',linestyle='-')
    plot(tz,sy,label='Sy[m]',color='red' ,linestyle='-')
    vscale=axis()[3]*0.1
    viseoy = [x*vscale for x in vis_abszisse]
    for i in range(stop_viseo,len(vis_ordinate)): viseoy[i] = 0.   # stop lattice plotting
    plot(vis_ordinate,viseoy,label='',color='black')
    plot(vis_ordinate,vzero,color='black')
    legend(loc='lower right',fontsize='x-small')

    show()

def display1(functions):
    """
    Plotting with longitudinal motion
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
    fighdr = 'lattice version = {}, input file = {}'.format(PARAMS['lattice_version'],PARAMS['input_file'])
    figure(fighdr,figsize=(width,height),facecolor='#eaecef',tight_layout=True)
    #-------------------- transverse X
    splot=subplot(311)
    splot.set_title('transverse x')
    plot(z,bx ,label=r'$\sigma$ [m]',color='green')
    plot(tz,cx,label='C [m]',color='blue',linestyle=':')
    plot(tz,sx,label='S [m]',color='red' ,linestyle=':')
    vscale=axis()[3]*0.4
    viseox = [x*vscale for x in vis_abszisse]
    for i in range(stop_viseo,len(vis_ordinate)): viseox[i] = 0.   # stop lattice plotting
    plot(vis_ordinate,viseox,label='',color='black')
    plot(vis_ordinate,vzero,color='black')
    legend(loc='lower right',fontsize='x-small')
    #-------------------- transverse Y
    splot=subplot(312)
    splot.set_title('transverse y')
    plot(z,by ,label=r'$\sigma$ [m]',color='green')
    plot(tz,cy,label='C [m]',color='blue',linestyle=':')
    plot(tz,sy,label='S [m]',color='red' ,linestyle=':')
    vscale=axis()[3]*0.4
    viseoy = [x*vscale for x in vis_abszisse]
    for i in range(stop_viseo,len(vis_ordinate)): viseoy[i] = 0.   # stop lattice plotting
    plot(vis_ordinate,viseoy,label='',color='black')
    plot(vis_ordinate,vzero,color='black')
    legend(loc='lower right',fontsize='x-small')
    #-------------------- longitudinal dPhi, dW/W
    ax_l=subplot(313)
    ax_l.set_title('longitudinal')
    ax_l.set_ylabel(r"$\Delta\phi$ [deg]")
    ax_l.tick_params(axis='y', colors='green')
    ax_l.yaxis.label.set_color('green')
    ax_l.plot(tz,cz,color='green',linestyle=':')
    ax_l.plot(tz,sz,color='green')
    vscale=ax_l.axis()[3]*0.7
    viseoz = [x*vscale for x in vis_abszisse]
    for i in range(stop_viseo,len(vis_ordinate)): viseoz[i] = 0.   # stop lattice plotting
    ax_l.plot(vis_ordinate[0:stop_viseo],viseoz[0:stop_viseo],label='',color='black')
    ax_l.plot(vis_ordinate,vzero,color='green',linestyle='--')

    ax_r = ax_l.twinx()
    ax_r.set_ylabel(r'$\Delta$w/w [%]')
    ax_r.tick_params(axis='y', colors='red')
    ax_r.yaxis.label.set_color('red')
    ax_r.plot(tz,cdw,color='red',linestyle=':')
    ax_r.plot(tz,sdw,color='red')
    ax_r.plot(vis_ordinate,vzero,color='red', linestyle='--')

    show()

# -------------------------START here
def simulation(filepath):
    lattice = parse_yaml_and_fabric(filepath)
    # Energie Konfiguration hier (SUPER WICHTIG)
    soll_track = track_soll(lattice)
    # sys.exit(9)
    lattice.stats(soll_track)          ## count elements and other statistics
    # ganzer Beschleuniger, Anfangswerte, Summary, etc...
    lattice.cell(closed=FLAGS['periodic'])
    collect_data_for_summary(lattice)    ## summary

    # Grafik & Lösungen
    lat_plot = lattice.lattice_plot_function()
    KVprint_flag = FLAGS['KVprint']
    if not KVprint_flag: print('CALCULATE C+S TRAJECTORIES')
    resolution = 23
    (c_like,s_like) = lattice.cs_traj(steps=resolution)       # calc sin- and cos-like trajectories
    if FLAGS['sigma']:
        if not KVprint_flag: print('CALCULATE SIGMA')
        sigma = lattice.sigma_functions(steps=resolution)     # calc. beamsize from sigma-matrix
    elif not FLAGS['sigma']:
        if not KVprint_flag: print('CALCULATE TWISS')
        twiss = lattice.twiss_functions(steps=resolution)     # calc. beamsize from beta-matrix
        sigma = [(x[0],sqrt(x[1]*PARAMS['emitx_i']),sqrt(x[2]*PARAMS['emity_i'])) for x in twiss]
    if KVprint_flag:
        dictprnt(PARAMS,text='PARAMS',njust=1)
    else:
        dictprnt(SUMMARY,text='summary')         # summary
        display((sigma,c_like,s_like,lat_plot))  # plot!!!

if __name__ == '__main__':
    filepath = 'yml/ref_run.yml'     ## the default input file (YAML syntax)
    filepath = 'yml/work.yml'

    if len(sys.argv) == 2:
        filepath = sys.argv[1]
    simulation(filepath)
