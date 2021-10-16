##!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
___version___='v8.0.6a2'
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
#todo: Calulate cell phase advance sigma by integration
#todo: handle exceptions speziel ValueError   - more or less done
#todo: use normalized emittances ?
#todo: waccept results are global - each node should carry its own
#todo: make new simu_manual.tex, README.md, check conversions.tex
#todo: rework the KVout - done in parts
#todo: rework verbose printing levels
#todo: C.K.Allen's matrices which are XAL as well?
#todo: slices as sub-lattice attribute to thick element

import sys
import os
# import subprocess
# from math import sqrt

# for PyQt
# import PyQt5             # works on native W10 but not on WSL2 as docker container
# import matplotlib
# matplotlib.use("Qt5Agg") # works on native W10 but not on WSL2 as docker container

# for Tk
import tkinter             # works on native W10
import matplotlib
matplotlib.use("TkAgg")    # works on native W10

import matplotlib.pyplot as plt
# plt.ioff()    # interactive mode off need this?
# plt.on()      # interactive mode on need this ?
# from matplotlib.patches import Ellipse

from setutil import PARAMS,FLAGS,SUMMARY,dictprnt,DEBUG_OFF,DEBUG_ON
from setutil import collect_data_for_summary, waccept, elli_sxy_action
from lattice_generator import factory
from tracker import track_soll

import bucket_size

DEBUG_MODULE   = DEBUG_OFF
DEBUG_LATTICE  = DEBUG_OFF

def bucket(*args):
    bucket_size.bucket()

def display0(*args):
    """
    C&S-Tracks w/o longitudinal motion
    """
    #----------*----------*   # unpack
    sigma_fun = args[0]
    cos_like  = args[1]
    sin_like  = args[2]
    lat_plot  = args[3]
    #-------------------- Bahnkoordinate (z)
    z     = [sigma_fun(i,'s')      for i in range(sigma_fun.nbpoints)]
    sgx   = [sigma_fun(i,'sigmax') for i in range(sigma_fun.nbpoints)]
    sgy   = [sigma_fun(i,'sigmay') for i in range(sigma_fun.nbpoints)]
#    zero  = [0.                    for i in range(sigma_fun.nbpoints)]
    #-------------------- trajectories (tz)
    tz=  [cos_like(i,'s')     for i in range(cos_like.nbpoints)]
    cx=  [cos_like(i,'cx')    for i in range(cos_like.nbpoints)]
#    cxp= [cos_like(i,'cxp')   for i in range(cos_like.nbpoints)]
    cy=  [cos_like(i,'cy')    for i in range(cos_like.nbpoints)]
#    cyp= [cos_like(i,'cyp')   for i in range(cos_like.nbpoints)]
#    cz=  [cos_like(i,'cz')    for i in range(cos_like.nbpoints)]
#    cdp= [cos_like(i,'cdp')   for i in range(cos_like.nbpoints)]

    sx=  [sin_like(i,'sx')    for i in range(sin_like.nbpoints)]
#    sxp= [sin_like(i,'sxp')   for i in range(sin_like.nbpoints)]
    sy=  [sin_like(i,'sy')    for i in range(sin_like.nbpoints)]
#    syp= [sin_like(i,'syp')   for i in range(sin_like.nbpoints)]
#    sz=  [sin_like(i,'sz')    for i in range(sin_like.nbpoints)]
#    sdp= [sin_like(i,'sdp')   for i in range(sin_like.nbpoints)]
    #-------------------- lattice viseo
    stop_viseox  = 5                  # stop viseo plot after so many [m]
    stop_viseoy  = 5                  # stop viseo plot after so many [m]
    vzero        = [0.                  for i in range(lat_plot.nbpoints)] # zero line
    vis_abszisse = [lat_plot(i,'s')     for i in range(lat_plot.nbpoints)]
    vis_ordinate = [lat_plot(i,'viseo') for i in range(lat_plot.nbpoints)]
    #-------------------- figure frame
    width=14; height=7.6
    plt.figure(num=0,figsize=(width,height),facecolor='#eaecef',tight_layout=False)

    #-------------------- transverse X
    splot211=plt.subplot(211)
    splot211.set_title('transverse x')
    plt.plot(z,sgx ,label=r'$\sigma$ [m]',color='green')
    plt.plot(tz,cx ,label='Cx[m]', color='blue',linestyle='-')
    # plt.plot(tz,cxp,label="Cx'[m]",color='blue',linestyle=':')
    plt.plot(tz,sx, label='Sx[m]', color='red' ,linestyle='-')
    # plt.plot(tz,sxp,label="Sx'[m]",color='red' ,linestyle=':')
    vscale=plt.axis()[3]*0.1
    viseox = [x*vscale for x in vis_ordinate]
    for i,s in enumerate(vis_abszisse):
        if s > stop_viseox:
            viseox[i] = 0.
    plt.plot(vis_abszisse,viseox,label='',color='black')
    plt.plot(vis_abszisse,vzero,color='black')
    plt.legend(loc='lower right',fontsize='x-small')

    #-------------------- transverse Y
    splot212=plt.subplot(212)
    splot212.set_title('transverse y')
    plt.plot(z,sgy ,label=r'$\sigma$ [m]',color='green')
    plt.plot(tz,cy, label='Cy[m]', color='blue',linestyle='-')
    # plt.plot(tz,cyp,label="Cy'[m]",color='blue',linestyle=':')
    plt.plot(tz,sy, label='Sy[m]', color='red' ,linestyle='-')
    # plt.plot(tz,syp,label="Sy'[m]",color='red' ,linestyle=':')
    vscale=plt.axis()[3]*0.1
    viseoy = [x*vscale for x in vis_ordinate]
    for i,s in enumerate(vis_abszisse):
        if s > stop_viseoy:
            viseoy[i] = 0.
    plt.plot(vis_abszisse,viseoy,label='',color='black')
    plt.plot(vis_abszisse,vzero,color='black')
    plt.legend(loc='lower right',fontsize='x-small')

def display1(*args):
    """
    C&S-Tracks with longitudinal motion
    """
    #-------------------- unpack
    sigma_fun = args[0]
    cos_like  = args[1]
    sin_like  = args[2]
    lat_plot  = args[3]
    ape_plot  = args[4]
    #-------------------- sigma functions
#    zero  = [0.                    for i in range(sigma_fun.nbpoints)] # zero line
    z     = [sigma_fun(i,'s')      for i in range(sigma_fun.nbpoints)] # Abszisse
    sgx   = [sigma_fun(i,'sigmax')*1.e3 for i in range(sigma_fun.nbpoints)] # envelope (sigma-x)
    sgy   = [sigma_fun(i,'sigmay')*1.e3 for i in range(sigma_fun.nbpoints)] # envelope (sigma-y)
    #-------------------- trajectories
    z1=  [cos_like(i,'s')          for i in range(cos_like.nbpoints)]
    cx=  [cos_like(i,'cx')*1.e3    for i in range(cos_like.nbpoints)]
#    cxp= [cos_like(i,'cxp')*1.e3   for i in range(cos_like.nbpoints)]
    cy=  [cos_like(i,'cy')*1.e3    for i in range(cos_like.nbpoints)]
#    cyp= [cos_like(i,'cyp')*1.e3   for i in range(cos_like.nbpoints)]
    cz=  [cos_like(i,'cz')         for i in range(cos_like.nbpoints)]
    cdp= [cos_like(i,'cdp')        for i in range(cos_like.nbpoints)]

    z2=  [sin_like(i,'s')          for i in range(sin_like.nbpoints)]
    sx=  [sin_like(i,'sx')*1.e3    for i in range(sin_like.nbpoints)]
#    sxp= [sin_like(i,'sxp')*1.e3   for i in range(sin_like.nbpoints)]
    sy=  [sin_like(i,'sy')*1.e3    for i in range(sin_like.nbpoints)]
#    syp= [sin_like(i,'syp')*1.e3   for i in range(sin_like.nbpoints)]
    sz=  [sin_like(i,'sz')         for i in range(sin_like.nbpoints)]
    sdp= [sin_like(i,'sdp')        for i in range(sin_like.nbpoints)]
    #-------------------- lattice viseo
    stop_viseox = 5                  # stop viseo plot after so many [m]
    stop_viseoy = 5                  # stop viseo plot after so many [m]
    stop_viseoz = 5                  # stop viseo plot after so many [m]
    vzero        = [0.                           for i in range(lat_plot.nbpoints)] # zero line
    vis_abszisse = [lat_plot(i,'s')              for i in range(lat_plot.nbpoints)]
    vis_ordinate = [lat_plot(i,'viseo')          for i in range(lat_plot.nbpoints)]
    ape_abszisse = [ape_plot(i,'s')              for i in range(ape_plot.nbpoints)]
    ape_ordinate = [ape_plot(i,'aperture')*1.e3  for i in range(ape_plot.nbpoints)]
    #-------------------- figure frame
    width=14; height=7.6
    # fighdr = 'lattice version = {}, input file = {}'.format(PARAMS['lattice_version'],PARAMS['input_file'])
    plt.figure(num=1,figsize=(width,height),facecolor='#eaecef',tight_layout=False)

    #-------------------- transverse X tracks
    splot311=plt.subplot(311)
    splot311.set_title('transverse x')
    # mapping box
    splot311.text(0.01, 1.1, PARAMS['mapping'],transform=splot311.transAxes,fontsize=8,bbox=dict(boxstyle='round',facecolor='wheat',alpha=0.5),verticalalignment='top')
    plt.plot(z,sgx ,label=r'$\sigma$ [mm]',color='green')
    plt.plot(z1,cx, label="C  [mm]",color='blue',linestyle='-')
    # plt.plot(z1,cxp,label="C' [mr]",color='blue',linestyle=':')
    plt.plot(z2,sx, label="S  [mm]",color='red' ,linestyle='-')
    # plt.plot(z2,sxp,label="S' [mr]",color='red' ,linestyle=':')
    # lattice elements
    vscale=plt.axis()[3]*0.4
    viseox = [x*vscale for x in vis_ordinate]
    # stop lattice plotting after stop_viseo meters
    for i,s in enumerate(vis_abszisse):
        if s > stop_viseox:
            viseox[i] = 0.
    plt.plot(vis_abszisse,viseox,label='',color='black')
    plt.plot(vis_abszisse,vzero,color='black')
    # apertures
    if FLAGS['useaper']:
        plt.plot(ape_abszisse,ape_ordinate,linestyle='-.')
        N = PARAMS['nbsigma']
        sgx = [i*N for i in sgx]
        #label = F'{N:1}$\sigma$ [mm]'
        label = '{:1}$\sigma$ [mm]'.format(N)
        plt.plot(z,sgx ,label=label,color='green',linestyle=':')
    plt.legend(loc='lower right',fontsize='x-small')

    #-------------------- transverse Y tracks
    splot312=plt.subplot(312)
    splot312.set_title('transverse y')
    plt.plot(z,sgy ,label=r'$\sigma$ [mm]',color='green')
    plt.plot(z1,cy, label="C  [mm]",color='blue',linestyle='-')
    # plt.plot(z1,cyp,label="C' [mr]",color='blue',linestyle=':')
    plt.plot(z2,sy, label="S  [mm]",color='red' ,linestyle='-')
    # plt.plot(z2,syp,label="S' [mr]",color='red' ,linestyle=':')
    # lattice elements
    vscale=plt.axis()[3]*0.4
    viseoy = [x*vscale for x in vis_ordinate]
    # stop lattice plotting after stop_viseo meters
    for i,s in enumerate(vis_abszisse):
        if s > stop_viseoy:
            viseoy[i] = 0.
    plt.plot(vis_abszisse,viseoy,label='',color='black')
    plt.plot(vis_abszisse,vzero,color='black')
    # apertures
    if FLAGS['useaper']:
        plt.plot(ape_abszisse,ape_ordinate,linestyle='-.')
        N = PARAMS['nbsigma']
        sgy = [i*N for i in sgy]
        plt.plot(z,sgy ,label=label,color='green',linestyle=':')
    plt.legend(loc='lower right',fontsize='x-small')

    #-------------------- longitudinal tracks z, dP/P
    # ax_l = left abszisse
    ax_l=plt.subplot(313)
    ax_l.set_title('longitudinal')
    ax_l.set_ylabel(r"z [mm]")
    ax_l.tick_params(axis='y', colors='green')
    ax_l.yaxis.label.set_color('green')
    ax_l.plot(z1,cz,color='green')
    ax_l.plot(z2,sz,color='green',linestyle=':')
    # ax_r = right abszisse
    ax_r = ax_l.twinx()
    ax_r.set_ylabel(r'$\Delta$p/p [%]')
    ax_r.tick_params(axis='y', colors='red')
    ax_r.yaxis.label.set_color('red')
    ax_r.plot(z2,cdp,color='red')
    ax_r.plot(z2,sdp,color='red',linestyle=':')
    ax_r.plot(vis_abszisse,vzero,color='red', linestyle='--')
    # lattice elements
    vscale=ax_l.axis()[3]*0.7
    viseoz = [x*vscale for x in vis_ordinate]
    # stop lattice plotting after stop_viseo meters
    for i,s in enumerate(vis_abszisse):
        if s > stop_viseoz:
            viseoz[i] = 0.
    ax_l.plot(vis_abszisse,viseoz,label='',color='black')
    ax_l.plot(vis_abszisse,vzero,color='green',linestyle='--')

def display2(*args):
    elli_sxy_action(on_injection=True)

#                            |----------------------- |
# -------------------------  | everything starts here |
#                            |----------------------- |
def simulation(filepath):
    def display(*functions):
        plots   = []
        if FLAGS['csTrak'] and FLAGS['dWf'] == 0:
            plots.append(display0) # C&S tracks {x,y}
        elif FLAGS['csTrak'] and FLAGS['dWf'] == 1:
            plots.append(display1) # C&S tracks {x,y,z}
            if FLAGS['bucket']:
                plots.append(bucket) # separatrix
        if FLAGS['pspace']:
            plots.append(display2)   # pspace

        # standard plots
        if len(plots) != 0:
            print('PREPARE DISPLAY')
            [plot(*functions) for plot in plots]
    #todo: what about markers plots ?
        # lattice.marker_actions()
        # next not needed for jupyter
        plt.show()

    #----------------------------------------------
    # STEP 1: parse input file and create a lattice
    #----------------------------------------------
    lattice = factory(filepath)

    if DEBUG_LATTICE(): lattice.show_linkage()

    #----------------------------------------------
    # STEP 2: configure elements for energy increase
    #----------------------------------------------
    soll_track = track_soll(lattice)

    if DEBUG_LATTICE(): lattice.show_linkage()

    #print(F'FINAL kinetic energy {lattice.seq[-1].particle.T} [MeV]')
    print('FINAL kinetic energy {} [MeV]'.format(lattice.seq[-1].particle.T))

    #----------------------------------------------
    # STEP 3: calculate longitudinal paramters at entrance
    #----------------------------------------------
    waccept(lattice.first_gap)

    #----------------------------------------------
    # STEP 4: count elements and make other statistics
    #----------------------------------------------
    lattice.stats(soll_track)

    #----------------------------------------------
    # STEP 5: beam dynamics full accelerator: initial values, etc...
    #----------------------------------------------
    lattice.cell(closed = FLAGS['periodic'])

    #----------------------------------------------
    # STEP 6:collect results and display them
    #----------------------------------------------
    collect_data_for_summary(lattice)

    # results
    kv_only = FLAGS['KVout']
    if kv_only:
        kv = {}
        for key in PARAMS:
            kv[key] = PARAMS[key]
        for key in SUMMARY:
            kv[key] = SUMMARY[key]
        dictprnt(kv,text='KV',njust=1)
    else:
        steps = 10
        # show summary
        dictprnt(SUMMARY,text='summary')
        # generate lattice plot
        (lat_plot, ape_plot) = lattice.lattice_plot_functions()
        # track sin- and cos-like trajectories
        (c_like,s_like) = lattice.cs_traj(steps = steps)
        # calculate envelope functions
        sigma_fun = lattice.sigmas(steps = steps)
        # make plots of functionsa
        display(sigma_fun,c_like,s_like,lat_plot,ape_plot)

if __name__ == '__main__':
    print('simu.py {} on python {}.{}.{} on {}'.format(___version___,sys.version_info.major,sys.version_info.minor,sys.version_info.micro,sys.platform))

    macros_file = ''
    template_file = ''
    # preset files for launch with  m4
    run_version   = '20.02.2019'
    input_file    = 'yml/simuIN.yml'       # default input file        (UNIX EOL=LF)


    if len(sys.argv) == 2:
        input_file    = sys.argv[1]
    else:
        if sys.platform   == 'win32':
            # launch .bat script
            command = 'yml\m4_simu.bat '+run_version
            print(command)
        elif sys.platform == 'darwin' or sys.platform.startswith('linux'):
            macros_file   = 'yml/macros_'+run_version+'.sh'
            template_file = 'yml/tmpl_'+run_version+'.yml'
            # launch bash
            command = 'chmod +x {}'.format(macros_file)
            command = "{0};{1} {2} {3}".format(command,macros_file,template_file, input_file)
        else:
            print('wrong platform')
            sys.exit(1)
        print('run Version {0}\n   macros=macros_{0}\n   template=tmpl_{0}\n   input={1}'.format(run_version,input_file))
        # print(command)
        os.system(command)

    # start the run
    simulation(input_file)

