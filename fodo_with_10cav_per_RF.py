#!/Users/klotz/SIMULINAC_env/bin/python
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
from setup import CONF,SUMMARY,dictprnt
from pylab import plot,show,legend,figure,subplot,axis
from math import sqrt
from fileLoader import read_yaml_and_parse

def display(functions):
    if CONF['dWf'] == 0:
        display0(functions)
    else:
        display1(functions)

def display0(functions):          ## plotting w/o longitudinal motion
    #----------*----------*   # unpack
    beta_fun = functions[0]
    cos_like = functions[1]
    sin_like = functions[2]
    emix=CONF['emitx_i']  # emittance @ entrance
    emiy=CONF['emity_i']  # emittance @ entrance
    #----------*----------*   # bahnkoordinate z
    z   = [ x[0] for x in beta_fun]
    #----------*----------*
    bx  = [ sqrt(x[1]*emix) for x in beta_fun]    # envelope (beta-x)
    by  = [ sqrt(x[2]*emiy) for x in beta_fun]    # envelope (beta-y)
#     bxn = [-x for x in bx]    # beta-x (negatif)
#     byn = [-x for x in by]    # beta-y (negatif)
    #----------*----------*   # trajectories
    cx = [x[0] for x in cos_like]   # cos-like-x
    cy = [x[2] for x in cos_like]   # cos-like-y
    cz = [x[4] for x in cos_like]   # cos-like-z
    cdw= [x[5] for x in cos_like]   # cos-like-dw/w
    sx = [x[0] for x in sin_like]   # sin-like-x
    sy = [x[2] for x in sin_like]   # sin-like-x
    sz = [x[4] for x in sin_like]   # sin-like-z
    sdw= [x[5] for x in sin_like]   # sin-like-dw/w
    #----------*----------*   # figure frame
    viseo = [x[3] for x in beta_fun]
    zero  = [0.   for x in beta_fun]# zero line
    width=20; height=12
    figure(SUMMARY['lattice_version'],figsize=(width,height))
    # figure(SUMMARY['lattice_version'])
    #----------*----------*   # transverse X
    splot=subplot(211)
    splot.set_title('transverse x')
    plot(z,bx ,label=r'$\sigma$ [m]',color='green')
#     plot(z,bxn,label='',color='green')
    plot(z,cx,label='Cx[m]',color='blue',linestyle='-.')
    plot(z,sx,label='Sx[m]',color='red' ,linestyle='-.')
    vscale=axis()[3]*0.1
    viseox = [x*vscale for x in viseo]
    plot(z,viseox,label='',color='black')
    plot(z,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    #----------*----------*   # transverse Y
    splot=subplot(212)
    splot.set_title('transverse y')
    plot(z,by ,label=r'$\sigma$ [m]',color='green')
#     plot(z,byn,label='',color='green')
    plot(z,cy,label='Cy[m]',color='blue',linestyle='-.')
    plot(z,sy,label='Sy[m]',color='red' ,linestyle='-.')
    vscale=axis()[3]*0.1
    viseoy = [x*vscale for x in viseo]
    plot(z,viseoy,label='',color='black')
    plot(z,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    #----------*----------*
    show(block=True)

def display1(functions):          ## plotting with longitudinal motion
    #----------*----------*   # unpack
    beta_fun = functions[0]
    cos_like = functions[1]
    sin_like = functions[2]
    emix=CONF['emitx_i']  # emittance @ entrance
    emiy=CONF['emity_i']  # emittance @ entrance
    #----------*----------*   # bahnkoordinate z
    z   = [ x[0] for x in beta_fun]
    #----------*----------*
    bx  = [ sqrt(x[1]*emix) for x in beta_fun]    # envelope (beta-x)
    by  = [ sqrt(x[2]*emiy) for x in beta_fun]    # envelope (beta-y)
#     bxn = [-x for x in bx]    # beta-x (negatif)
#     byn = [-x for x in by]    # beta-y (negatif)
    #----------*----------*   # trajectories
    cx = [x[0] for x in cos_like]   # cos-like-x
    cy = [x[2] for x in cos_like]   # cos-like-y
    cz = [x[4] for x in cos_like]   # cos-like-z
    cdw= [x[5] for x in cos_like]   # cos-like-dw/w
    sx = [x[0] for x in sin_like]   # sin-like-x
    sy = [x[2] for x in sin_like]   # sin-like-x
    sz = [x[4] for x in sin_like]   # sin-like-z
    sdw= [x[5] for x in sin_like]   # sin-like-dw/w
    #----------*----------*   # figure frame
    viseo = [x[3] for x in beta_fun]
    zero  = [0.   for x in beta_fun]# zero line
    width=20; height=12
    figure(SUMMARY['lattice_version'],figsize=(width,height))
    # figure(SUMMARY['lattice_version'])
    #----------*----------*   # transverse X
    splot=subplot(311)
    splot.set_title('transverse x')
    plot(z,bx ,label=r'$\sigma$ [m]',color='green')
#     plot(z,bxn,label='',color='green')
    plot(z,cx,label='Cx[m]',color='blue',linestyle='-.')
    plot(z,sx,label='Sx[m]',color='red' ,linestyle='-.')
    vscale=axis()[3]*0.1
    viseox = [x*vscale for x in viseo]
    plot(z,viseox,label='',color='black')
    plot(z,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    #----------*----------*   # transverse Y
    splot=subplot(312)
    splot.set_title('transverse y')
    plot(z,by ,label=r'$\sigma$ [m]',color='green')
#     plot(z,byn,label='',color='green')
    plot(z,cy,label='Cy[m]',color='blue',linestyle='-.')
    plot(z,sy,label='Sy[m]',color='red' ,linestyle='-.')
    vscale=axis()[3]*0.1
    viseoy = [x*vscale for x in viseo]
    plot(z,viseoy,label='',color='black')
    plot(z,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    #----------*----------*   # longitudinal dPhi, dW/W
    ax_l=subplot(313)
    ax_l.set_title('longitudinal z')
    ax_l.set_ylabel(r"$\Delta\phi$ [deg]")
    ax_l.tick_params(axis='y', colors='green')
    ax_l.yaxis.label.set_color('green')
    ax_l.plot(z,cz,label=r"$\Delta\phi$"  ,color='green')
    ax_l.plot(z,sz,color='green')
    vscale=ax_l.axis()[3]*0.1
    viseoz = [x*vscale for x in viseo]
    ax_l.plot(z,viseoz,label='',color='black')
    ax_l.plot(z,zero,color='black')

    ax_r = ax_l.twinx()
    ax_r.set_ylabel(r'$\Delta$w/w [%]')
    ax_r.tick_params(axis='y', colors='red')
    ax_r.yaxis.label.set_color('red')
    ax_r.plot(z,cdw,label=r'$\Delta$w/w',color='red')
    ax_r.plot(z,sdw,color='red')
    ax_r.plot(z,zero,color='red', linestyle='--')
    #----------*----------*
    show(block=True)

def loesung(filepath):                   ## total classic FODO lattice (1st result, used as reference!)
    super_cell = read_yaml_and_parse(filepath)
    #-----------------------------------------
    # Rechne: ganze Zelle und Anfangswerte
    mcell,betax,betay = super_cell.cell(closed=CONF['periodic'])
    dictprnt(SUMMARY,text='summary')
    #-----------------------------------------
    # Zeige: Grafik Lösungen als Funktion von (s)
    functions = super_cell.functions(30)
    display(functions)

if __name__ == '__main__':
    import sys, os
    directory = os.path.dirname(__file__)
#     filepath = directory+'/fodo_with_10cav_per_RF.yml'       ## the default input file (YAML syntax)
    filepath = directory+'/fodo_with_10cav_per_RF(1).yml'       ## the default input file (YAML syntax)
    if len(sys.argv) == 2:
        filepath = sys.argv[1]
    loesung(filepath)
