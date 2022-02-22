##!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
___version___='v10.0.1'
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
#TODO: Calulate cell phase advance sigma by integration
#TODO: handle exceptions speziel ValueError   - more or less done
#TODO: use normalized emittances ?
#TODO: waccept results are global - each node should carry its own
#TODO: make new simu_manual.tex, README.md, check conversions.tex
#TODO: rework the KVout - done in parts
#TODO: rework verbose printing levels
#TODO: C.K.Allen's matrices which are XAL as well?
#TODO: slices as sub-lattice attribute to thick element
#TODO: for tracker: Plot a confidence ellipse of a two-dimensional dataset: https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html#sphx-glr-gallery-statistics-confidence-ellipse-py
#TODO: Covariance Ellipse see https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html
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
# from matplotlib.patches import Ellipse
import pprint, inspect

import bucket_size
from setutil import PARAMS,FLAGS,SUMMARY,dictprnt,waccept
from setutil import collect_data_for_summary, show_data_from_elements
from lattice_generator import factory
from PsMarkerAgent import ellipse_plot
# from tracker import track_soll
from pargs import pargs
from lattice_parser2 import parse as getParseResult

def PRINT_PRETTY(obj=None):
    file = inspect.stack()[0].filename
    print(F'DEBUG_ON[{file}] ==> ',end="")
    if obj != None: pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
    return True
def PASS(obj=None):
    return False
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

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
    # vscale=plt.axis()[3]*0.1
    # viseox = [x*vscale for x in vis_ordinate]
    # for i,s in enumerate(vis_abszisse):
    #     if s > stop_viseox:
    #         viseox[i] = 0.
    # plt.plot(vis_abszisse,viseox,label='',color='black')
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
    # for i,s in enumerate(vis_abszisse):
    #     if s > stop_viseoy:
    #         viseoy[i] = 0.
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
    vzero        = [0.                           for i in range(lat_plot.nbpoints)] # zero line
    vis_abszisse = [lat_plot(i,'s')              for i in range(lat_plot.nbpoints)]
    vis_ordinate = [lat_plot(i,'viseo')          for i in range(lat_plot.nbpoints)]
    ape_abszisse = [ape_plot(i,'s')              for i in range(ape_plot.nbpoints)]
    ape_ordinate = [ape_plot(i,'aperture')*1.e3  for i in range(ape_plot.nbpoints)]
    #-------------------- figure frame
    width=14; height=7.6
    # fighdr = 'lattice version = {}, input file = {}'.format(PARAMS['lattice_version'],PARAMS['input_file'])
    fig = plt.figure(num=1,figsize=(width,height),facecolor='#eaecef',tight_layout=False)

    #-------------------- transverse X tracks
    splot311=plt.subplot(311)
    # splot311=plt.subplot(10,1,(1,3))
    splot311.set_title('transverse x')
    # mapping box
    splot311.text(0.01, 1.1, PARAMS['mapping'],transform=splot311.transAxes,fontsize=8,bbox=dict(boxstyle='round',facecolor='wheat',alpha=0.5),verticalalignment='top')
    plt.plot(z,sgx ,label=r'$\sigma$ [mm]',color='green')
    plt.plot(z1,cx, label="C  [mm]",color='blue',linestyle='-')
    # plt.plot(z1,cxp,label="C' [mr]",color='blue',linestyle=':')
    plt.plot(z2,sx, label="S  [mm]",color='red' ,linestyle='-')
    # plt.plot(z2,sxp,label="S' [mr]",color='red' ,linestyle=':')
    vscale=splot311.axis()[3]*0.25
    viseoz = [x*vscale for x in vis_ordinate]
    plt.plot(vis_abszisse,viseoz,label='',color='black')
    plt.plot(vis_abszisse,vzero,color='green',linestyle='--')
    # apertures
    if FLAGS['useaper']:
        plt.plot(ape_abszisse,ape_ordinate,linestyle='-.')
        N = PARAMS['nbsigma']
        sgx = [i*N for i in sgx]
        #label = F'{N:1}$\sigma$ [mm]'
        label = '{:1}$\sigma$ [mm]'.format(N)
        plt.plot(z,sgx ,label=label,color='green',linestyle=':')
    # zero line
    splot311.plot(vis_abszisse,vzero,color='green',linestyle='--')
    plt.legend(loc='lower right',fontsize='x-small')

    #-------------------- transverse Y tracks
    splot312=plt.subplot(312)
    # splot312=plt.subplot(10,1,(4,6))
    splot312.set_title('transverse y')
    plt.plot(z,sgy ,label=r'$\sigma$ [mm]',color='green')
    plt.plot(z1,cy, label="C  [mm]",color='blue',linestyle='-')
    # plt.plot(z1,cyp,label="C' [mr]",color='blue',linestyle=':')
    plt.plot(z2,sy, label="S  [mm]",color='red' ,linestyle='-')
    vscale=splot312.axis()[3]*0.25
    viseoz = [x*vscale for x in vis_ordinate]
    plt.plot(vis_abszisse,viseoz,label='',color='black')
    plt.plot(vis_abszisse,vzero,color='green',linestyle='--')
    # apertures
    if FLAGS['useaper']:
        plt.plot(ape_abszisse,ape_ordinate,linestyle='-.')
        N = PARAMS['nbsigma']
        sgy = [i*N for i in sgy]
        plt.plot(z,sgy ,label=label,color='green',linestyle=':')
    # zero line
    splot312.plot(vis_abszisse,vzero,color='green',linestyle='--')
    plt.legend(loc='lower right',fontsize='x-small')

    #-------------------- longitudinal tracks z, dP/P
    # ax_l = left abszisse
    ax_l=plt.subplot(313)
    # ax_l=plt.subplot(10,1,(7,9))
    ax_l.set_title('longitudinal')
    ax_l.set_ylabel(r"z [mm]")
    ax_l.tick_params(axis='y', colors='green')
    ax_l.yaxis.label.set_color('green')
    ax_l.plot(z1,cz,label='C',color='green')
    ax_l.plot(z2,sz,label='S',color='green',linestyle=':')
    plt.legend(loc='lower left',fontsize='x-small')
    # ax_r = right abszisse
    ax_r = ax_l.twinx()
    ax_r.set_ylabel(r'$\Delta$p/p [%]')
    ax_r.tick_params(axis='y', colors='red')
    ax_r.yaxis.label.set_color('red')
    ax_r.plot(z2,cdp,label='C',color='red')
    ax_r.plot(z2,sdp,label='S',color='red',linestyle=':')
    ax_r.plot(vis_abszisse,vzero,color='red', linestyle='--')
    plt.legend(loc='lower right',fontsize='x-small')
    # lattice elements
    vscale=ax_l.axis()[3]*0.25
    viseoz = [x*vscale for x in vis_ordinate]
    ax_l.plot(vis_abszisse,viseoz,label='',color='black')
    ax_l.plot(vis_abszisse,vzero,color='green',linestyle='--')
def display2(*args):
    ellipse_plot(None,on_injection=True)
def lattice_check(lattice):
    for x in lattice.seq:
        DEBUG_ON(x.label)
        if x.type == 'QFth' : print(x.type,x.matrix)
        if x.type == 'QDth' : print(x.type,x.matrix)
def link_check(lattice):
    DEBUG_ON()
    lattice.show_linkage()
# ------- everything starts here ------- everything starts here ------- everything starts here ------- everything starts here
# ------- everything starts here ------- everything starts here ------- everything starts here ------- everything starts here
# ------- everything starts here ------- everything starts here ------- everything starts here ------- everything starts here
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
    #----------------------------------------------
    # STEP 1: parse input file and create a lattice
    #         with linked and energy adjusted nodes.
    #----------------------------------------------
    lattice = factory(filepath)
    if 0: lattice_check(lattice)
    if 0: link_check(lattice)
    descriptor = getParseResult().DESCRIPTOR  # get DESCRIPTOR from parsed results
    if descriptor != None: print(descriptor)
    #----------------------------------------------
    # STEP 2: configure elements for energy increase
    #----------------------------------------------
    # soll_track = track_soll(lattice, PARAMS['injection_energy'])
    if 0: lattice_check(lattice)
    if 0: link_check(lattice)
    print(F'FINAL kinetic energy {lattice.seq[-1].particle.tkin} [MeV]')   #TODO make property og latttice
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
    # STEP 6:collect results
    #----------------------------------------------
    collect_data_for_summary(lattice)
    #----------------------------------------------
    # STEP 7:display results and display them
    #----------------------------------------------
    kv_only = FLAGS['KVout']
    if kv_only:
        kv = {}
        for key in PARAMS:
            kv[key] = PARAMS[key]
        for key in SUMMARY:
            kv[key] = SUMMARY[key]
        dictprnt(kv,text='KV',njust=1)
    else:
        show_data_from_elements() #..................................show ELEMENT attributes
        dictprnt(SUMMARY,text='Summary') #...........................show summary
        (lat_plot, ape_plot) = lattice.lattice_plot_functions() #....generate lattice plot
        steps = 1
        (c_like,s_like) = lattice.cs_traj(steps=steps) #..............track sin- and cos-like trajectories
        sigma_fun       = lattice.sigmas(steps=steps) #.,,,.................calculate envelope functions
        display(sigma_fun,c_like,s_like,lat_plot,ape_plot) #.........make plots of functionsa
        if FLAGS['marker']: #........................................any MARKER with actions?
            # TODO need a function to filter/print lattice elements
            for node in lattice.seq:
                if node.type != "MRK": continue
                DEBUG_OFF(node.toString())
                node.do_actions()
        """ show all figures - (must be the only one!) """
        plt.show()
if __name__ == '__main__':
    print('simu.py {} on python {}.{}.{} on {}'.format(___version___,sys.version_info.major,sys.version_info.minor,sys.version_info.micro,sys.platform))

    # parse argv and normalize
    print("sys.argv: {}".format(sys.argv))
    Args = pargs(sys.argv)
    print('This run: input({}), template({}), macro({})'.format(Args['file'],Args['tmpl'],Args['macro']))

    input_file = Args['file']
    if sys.platform == 'win32':
        if Args['mode']   == 'no_m4':
            pass
        elif Args['mode'] == 'm4':
            command = 'yml\m4_launch.bat {} {} {}'.format(Args['file'],Args['tmpl'],Args['macro'])
            stat = os.system(command)
            if stat != 0:
                print('\nWARNING: system-command returned error - using standard "yml/simuIN.yml" without m4-preprocessing!')
                print('WARNING: system-command returned error - using standard "yml/simuIN.yml" without m4-preprocessing!')
                print('WARNING: system-command returned error - using standard "yml/simuIN.yml" without m4-preprocessing!\n')
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
        else:
            print('Internal error!')
            sys.exit(1)
    else:
        print('wrong platform')
        sys.exit(1)
    # run the simulation
    simulation(input_file)
    # simulation('yml/tmpl_25.10.2021_new.yml')
