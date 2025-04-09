#!python
# -*- coding: utf-8 -*-
__version__='v11.0.3'
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

#TODO: shift initial beam paramters to bunch.py or new beam.py module -has 1st priority
#TODO: Calulate cell phase advance sigma by integration - maybe
#TODO: use normalized emittances why?
#TODO: rework verbose printing levels - needed?
#TODO: C.K.Allen's matrices which are XAL as well? - don't know if better
#TODO: slices as sub-lattice attribute to thick element - too big a modification? advantage?
#TODO: Covariance Ellipse see https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html - done
#TODO: handle exceptions speziel ValueError - more or less done
#TODO: for tracker: plot confidence ellipse - used reference: https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html#sphx-glr-gallery-statistics-confidence-ellipse-py - done
"""
import matplotlib
# import PyQt5
# matplotlib.use("Qt5Agg")
# import tkinter
matplotlib.use("TkAgg")
import sys
import argparse
from setutil import DEBUG_ON,DEBUG_OFF,EKOO,wrapRED

import lattice_generator as LG
import math              as M
import matplotlib.pyplot as plt
import setutil           as UTIL
import elements          as ELM
import pandas            as pd
import bucket_size       as BKTSZ

def display0(*args):
    """ C&S-Tracks w/o longitudinal motion """
    #----------*----------*   # unpack
    twiss_func = args[0]
    cos_like   = args[1]
    sin_like   = args[2]
    lat_plot   = args[3]
    #-------------------- Bahnkoordinate (z)
    z     = [twiss_func(i,'s')    for i in range(twiss_func.nbpoints)]
    sgx   = [twiss_func(i,'sigx') for i in range(twiss_func.nbpoints)]
    sgy   = [twiss_func(i,'sigy') for i in range(twiss_func.nbpoints)]
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
    """ beta functions w/o longitudinal motion """
    #----------*----------*   # unpack
    twiss_func = args[0]
    lat_plot   = args[3]
    #-------------------- beta x,y & dispersion x
    s     = [twiss_func(i,'s')  for i in range(twiss_func.nbpoints)]
    bx    = [twiss_func(i,'bx') for i in range(twiss_func.nbpoints)]
    by    = [twiss_func(i,'by') for i in range(twiss_func.nbpoints)]
    dx    = [twiss_func(i,'dx') for i in range(twiss_func.nbpoints)]
    #-------------------- lattice viseo
    vzero        = [0.                  for i in range(lat_plot.nbpoints)] # zero line
    vis_abszisse = [lat_plot(i,'s')     for i in range(lat_plot.nbpoints)]
    vis_ordinate = [lat_plot(i,'viseo') for i in range(lat_plot.nbpoints)]
    #-------------------- figure frame
    width=14; height=7.6
    plt.figure(num=0,figsize=(width,height),facecolor='#eaecef',tight_layout=False)

    #-------------------- transverse X
    splot111=plt.subplot(111)
    splot111.set_title('beta functions')
    plt.plot(s,bx,  label=r'$\beta_x$ [m]', color='red',  linestyle='-') # beta x
    plt.plot(s,by,  label=r'$\beta_y$ [m]', color='blue', linestyle='-') # beta y
    plt.plot(s,dx,  label=r'$\eta_x$ [m]' , color='green',linestyle='-') # dispersion x
    vscale=splot111.axis()[3]*0.25
    viseoz = [x*vscale for x in vis_ordinate]
    plt.plot(vis_abszisse,viseoz,label='',color='black')
    plt.plot(vis_abszisse,vzero,color='black')
    plt.legend(loc='lower right',fontsize='x-small')
def display3(*args):
    """ C&S-Tracks with longitudinal motion """
    #-------------------- unpack
    twiss_fun = args[0]
    cos_like  = args[1]
    sin_like  = args[2]
    lat_plot  = args[3]
    ape_plot  = args[4]
    #-------------------- sigma functions
    #    zero  = [0.                    for i in range(sigma_fun.nbpoints)] # zero line
    z     = [twiss_fun(i,'s')         for i in range(twiss_fun.nbpoints)] # Abszisse
    sgx   = [twiss_fun(i,'sigx')*1.e3 for i in range(twiss_fun.nbpoints)] # envelope (sigma-x)
    sgy   = [twiss_fun(i,'sigy')*1.e3 for i in range(twiss_fun.nbpoints)] # envelope (sigma-y)
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
    splot311.text(0.01, 1.1,UTIL.FLAGS.get('mapping'),transform=splot311.transAxes,fontsize=8,bbox=dict(boxstyle='round',facecolor='wheat',alpha=0.5),verticalalignment='top')
    if UTIL.FLAGS['envelope']:
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
    if UTIL.FLAGS['useaper']:
        plt.plot(ape_abszisse,ape_ordinate,linestyle='-.')
        N = UTIL.PARAMS['nbsigma']
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
    if UTIL.FLAGS['envelope']:
        plt.plot(z,sgy ,label=r'$\sigma$ [mm]',color='green')
    plt.plot(z1,cy, label="C  [mm]",color='blue',linestyle='-')
    # plt.plot(z1,cyp,label="C' [mr]",color='blue',linestyle=':')
    plt.plot(z2,sy, label="S  [mm]",color='red' ,linestyle='-')
    vscale=splot312.axis()[3]*0.25
    viseoz = [x*vscale for x in vis_ordinate]
    plt.plot(vis_abszisse,viseoz,label='',color='black')
    plt.plot(vis_abszisse,vzero,color='green',linestyle='--')
    # apertures
    if UTIL.FLAGS['useaper']:
        plt.plot(ape_abszisse,ape_ordinate,linestyle='-.')
        N = UTIL.PARAMS['nbsigma']
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
    ax_r.plot(z1,cdp,label='C',color='red')
    ax_r.plot(z2,sdp,label='S',color='red',linestyle=':')
    ax_r.plot(vis_abszisse,vzero,color='red', linestyle='--')
    plt.legend(loc='lower right',fontsize='x-small')
    # lattice elements
    vscale=ax_l.axis()[3]*0.25
    viseoz = [x*vscale for x in vis_ordinate]
    ax_l.plot(vis_abszisse,viseoz,label='',color='black')
    ax_l.plot(vis_abszisse,vzero,color='green',linestyle='--')
def display4(*args):
    """ beta functions and synchrotron oscillations """
    #-------------------- unpack
    twiss_func = args[0]
    cos_like   = args[1]
    sin_like   = args[2]
    lat_plot   = args[3]
    #-------------------- beta x,y & dispersion x
    s     = [twiss_func(i,'s')    for i in range(twiss_func.nbpoints)] # Abszisse
    bx    = [twiss_func(i,'bx')   for i in range(twiss_func.nbpoints)] # beta x
    by    = [twiss_func(i,'by')   for i in range(twiss_func.nbpoints)] # beta y
    dx    = [twiss_func(i,'dx')   for i in range(twiss_func.nbpoints)] # dispersion x
    #-------------------- longitudinal trajectories
    z1=  [cos_like(i,'s')          for i in range(cos_like.nbpoints)]
    cz=  [cos_like(i,'cz')         for i in range(cos_like.nbpoints)]
    cdp= [cos_like(i,'cdp')        for i in range(cos_like.nbpoints)]

    z2=  [sin_like(i,'s')          for i in range(sin_like.nbpoints)]
    sz=  [sin_like(i,'sz')         for i in range(sin_like.nbpoints)]
    sdp= [sin_like(i,'sdp')        for i in range(sin_like.nbpoints)]
    #-------------------- lattice viseo
    vzero        = [0.                           for i in range(lat_plot.nbpoints)] # zero line
    vis_abszisse = [lat_plot(i,'s')              for i in range(lat_plot.nbpoints)]
    vis_ordinate = [lat_plot(i,'viseo')          for i in range(lat_plot.nbpoints)]
    #-------------------- figure frame
    width=14; height=7.6
    # fighdr = 'lattice version = {}, input file = {}'.format(PARAMS['lattice_version'],PARAMS['input_file'])
    fig = plt.figure(num=1,figsize=(width,height),facecolor='#eaecef',tight_layout=False)

    #-------------------- beta functions
    splot211=plt.subplot(211)
    splot211.set_title('beta x,y')
    # mapping box
    splot211.text(0.01, 1.1, UTIL.FLAGS.get('mapping'),transform=splot211.transAxes,fontsize=8,bbox=dict(boxstyle='round',facecolor='wheat',alpha=0.5),verticalalignment='top')
    # function plots
    plt.plot(s,bx,   label=r"$\beta$x  [m]",  color='black', linestyle='-')
    plt.plot(s,by,   label=r"$\beta$y  [m]",  color='red',   linestyle='-')
    plt.plot(s,dx,   label=r'$\eta_x$ [m]' ,  color='green', linestyle='-') # dispersion x
    vscale=splot211.axis()[3]*0.25
    viseoz = [x*vscale for x in vis_ordinate]
    plt.plot(vis_abszisse,viseoz,label='',color='black')
    plt.plot(vis_abszisse,vzero,color='green',linestyle='--')
    # zero line
    splot211.plot(vis_abszisse,vzero,color='green',linestyle='--')
    plt.legend(loc='lower right',fontsize='x-small')

    #-------------------- longitudinal tracks z, dP/P
    # ax_l = left abszisse
    ax_l=plt.subplot(212)
    # ax_l=plt.subplot(10,1,(7,9))
    ax_l.set_title('synchrotron oscillation')
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

def simulation(filepath):
    def display(*functions):
        # dispatch to different plots according to FLAG settings
        plots   = []
        if UTIL.FLAGS['dWf'] == 0:
            # accel OFF
            if UTIL.FLAGS['csTrak']:
                plots.append(display0)      # C&S tracks and sigmas {x,y}
            else:
                plots.append(display1)      # beta functions {x,y}
        else:
            # accel ON
            if UTIL.FLAGS['csTrak']:
                plots.append(display3)      # C&S tracks and sigmas {x,y}
            else:
                plots.append(display4)      # beta functions {x,y}

            if UTIL.FLAGS['bucket']:
                first_gap_node = lattice.first_gap
                if first_gap_node != None:
                    BKTSZ.bucket(first_gap_node)       # separatrix
                else:
                    print("No 1st rg-gap in lattice? Can't plot W-acceptance.")
        # make all plots
        if len(plots) != 0:
            print('PREPARE DISPLAY')
            [plot(*functions) for plot in plots]
    
    """ ------- everything starts here ------- everything starts here ------- everything starts here ------- everything starts here """
    #----------------------------------------------
    # STEP 1: parse input file and create a lattice
    #         with links and adjusted energy
    #         using ref_track in lattice.add_node.
    #         Set run-mode from FLAGS.
    #----------------------------------------------
    lattice = LG.factory(filepath)
    # descriptor
    UTIL.SUMMARY['Description'] = UTIL.PARAMS['descriptor']
    #----------------------------------------------
    # STEP 2: calculate run mode
    #----------------------------------------------
    UTIL.FLAGS['accON'] = lattice.accON
    # run-mode
    twoflag = (UTIL.FLAGS.get('accON'), UTIL.FLAGS.get('periodic'))
    if twoflag == (True,True):   mode= UTIL.RUN_MODE[0]
    if twoflag == (True,False):  mode= UTIL.RUN_MODE[1]
    if twoflag == (False,True):  mode= UTIL.RUN_MODE[2]
    if twoflag == (False,False): mode= UTIL.RUN_MODE[3]
    UTIL.FLAGS['mode'] = mode
    print(f'running in \'{UTIL.FLAGS["mode"]}\' mode')
    # print(wrapRED(f'\u26dd  FINAL kinetic energy {lattice.seq[-1].ref_track[EKOO]:.3f} [MeV] \u26dd'))
    # print(wrapRED(f'\u26dd  FINAL particle kinetic energy {lattice.last_gap.particle.tkin:.3f} [MeV] \u26dd'))
    print(wrapRED(f'\u26dd  FINAL particle kinetic energy {lattice.seq[-1].particle.tkin:.3f} [MeV] \u26dd'))
    #----------------------------------------------
    # STEP 3: count elements and make other statistics
    #----------------------------------------------
    lattice.make_label()
    lattice.make_matrix()
    stats = lattice.stats()
    UTIL.SUMMARY['nbof quadrupoles*']   = stats['quad_cntr']
    UTIL.SUMMARY['nbof cavities*']      = stats['cavity_cntr']
    UTIL.SUMMARY['Tkin (i,f) [MeV]']    = (stats['tki'],stats['tkf'])
    UTIL.SUMMARY['lattice length* [m]'] = stats['latt_length']
    #----------------------------------------------
    # STEP 4: beam dynamics full accelerator: initial values, etc...
    #----------------------------------------------
    res = lattice.cell(UTIL.FLAGS['periodic'])
    # Update PARAMS
    UTIL.PARAMS.update(res)
    #---------------------------------------------- 
    # STEP 5: lattice functions
    #----------------------------------------------
    steps = 10
    (c_like,s_like) = lattice.cs_traj(steps=steps)     #..........track sin- and cos-like trajectories
    twiss_func      = lattice.twiss_funcs(steps=steps) #.,,,......calculate envelope functions
    #----------------------------------------------
    # STEP 6: collect results
    #----------------------------------------------
    UTIL.collect_data_for_summary(lattice)
    #----------------------------------------------
    # STEP 7: ouput results
    #----------------------------------------------
    kv_only = UTIL.FLAGS['KVout']
    g_disp  = UTIL.FLAGS['GDisp']
    if kv_only:
        def flatten_dict(d):
            [flat_dict] = pd.json_normalize(d).to_dict(orient='records')
            return flat_dict
        def default():
            kv=dict(FLAGS=UTIL.FLAGS,PARAMS=UTIL.PARAMS,SUMMARY=UTIL.SUMMARY,ELEMENTS=UTIL.ELEMENTS)
            fkv = flatten_dict(kv)
            DEBUG_ON(fkv)
        def custom():
            kv={}
            kv['phase_advance'] = '{:.2f} {:.2f}'.format(M.degrees(UTIL.PHADVX),M.degrees(UTIL.PHADVY))
            kv['B\'f'] = UTIL.ELEMENTS['QF']["B'"]
            kv['B\'d'] = UTIL.ELEMENTS['QD']["B'"]
            print(kv)
        # out = default
        out = default
        out()
    else:
        UTIL.show_data_from_elements() #...................................show ELEMENT attributes
        if g_disp :
            (lat_plot, ape_plot) = lattice.lattice_plot_functions() #.....generate lattice plot
            display(twiss_func,c_like,s_like,lat_plot,ape_plot) #.........make plots of functions
        for node in lattice.seq: #....................................filter on Markers and invoke actions
            if not isinstance(node,ELM.MRK): continue
            DEBUG_OFF(node.toString())
            DEBUG_OFF(node.__dict__)
            node.do_action()
        UTIL.dictprnt(UTIL.SUMMARY,text='Summary') #............................show summary
        """ show all figures - (must be the only call to show!) """
        if g_disp :
            plt.show()

if __name__ == '__main__':
    # ArgumentParser puts result in 'args'
    parser = argparse.ArgumentParser()
    group  = parser.add_mutually_exclusive_group( )
    group.add_argument ("--file", default="simuIN.yml", help="lattice input-file")
    args = vars(parser.parse_args())
    print('simu.py {} on python {}.{}.{} on {}'.format(__version__,sys.version_info.major,sys.version_info.minor,sys.version_info.micro,sys.platform))

    # let's go. All  input is parsed...
    input_file = args['file']    # default simuIN.yml
    print('This run: input({}))'.format(input_file))

    if sys.platform == 'darwin' or sys.platform.startswith('linux') or sys.platform == 'win32':
        print(f'Platform {sys.platform}')
    else:
        print('wrong platform')
        sys.exit(1)
    # run the simulation
    simulation(input_file)
