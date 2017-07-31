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
from math import sqrt
from matplotlib.pyplot import plot,show,legend,figure,subplot,axis

from lattice_generator import parse_yaml_and_fabric
from bucket_size import bucket
from tracks import track_soll
from setutil import CONF,SUMMARY,Particle,Proton,dictprnt,collect_summary
from setutil import epsiz

def display(functions):
    if CONF['dWf'] == 0:
        display0(functions)
    else:
        display1(functions)
        # bucket()             # separatrix

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
    width=14; height=7.6
    figure(CONF['lattice_version'],figsize=(width,height),
                            facecolor='#eaecef',tight_layout=True)
    # figure(SUMMARY['lattice_version'])
    ## transverse X
    splot=subplot(211)
    splot.set_title('transverse x')
    plot(z,bx ,label=r'$\sigma$ [m]',color='green')
#     plot(z,bxn,label='',color='green')
    plot(z,cx,label='Cx[m]',color='blue',linestyle='-')
    plot(z,sx,label='Sx[m]',color='red' ,linestyle='-')
    vscale=axis()[3]*0.1
    viseox = [x*vscale for x in viseo]
    plot(z,viseox,label='',color='black')
    plot(z,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    ## transverse Y
    splot=subplot(212)
    splot.set_title('transverse y')
    plot(z,by ,label=r'$\sigma$ [m]',color='green')
#     plot(z,byn,label='',color='green')
    plot(z,cy,label='Cy[m]',color='blue',linestyle='-')
    plot(z,sy,label='Sy[m]',color='red' ,linestyle='-')
    vscale=axis()[3]*0.1
    viseoy = [x*vscale for x in viseo]
    plot(z,viseoy,label='',color='black')
    plot(z,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    #----------*----------*
    show()

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
    width=14; height=7.6
    figure(CONF['lattice_version'],figsize=(width,height),
                            facecolor='#eaecef',tight_layout=True)
    # figure(CONF['lattice_version'])
    ##  transverse X
    splot=subplot(311)
    splot.set_title('transverse x')
    plot(z,bx ,label=r'$\sigma$ [m]',color='green')
#     plot(z,bxn,label='',color='green')
    plot(z,cx,label='C [m]',color='blue',linestyle='-')
    plot(z,sx,label='S [m]',color='red' ,linestyle='-')
    vscale=axis()[3]*0.3
    viseox = [x*vscale for x in viseo]
    plot(z,viseox,label='',color='black')
    plot(z,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    ## transverse Y
    splot=subplot(312)
    splot.set_title('transverse y')
    plot(z,by ,label=r'$\sigma$ [m]',color='green')
#     plot(z,byn,label='',color='green')
    plot(z,cy,label='C [m]',color='blue',linestyle='-')
    plot(z,sy,label='S [m]',color='red' ,linestyle='-')
#     vscale=axis()[3]*0.1
    viseoy = [x*vscale for x in viseo]
    plot(z,viseoy,label='',color='black')
    plot(z,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    ## longitudinal dPhi, dW/W
    ax_l=subplot(313)
    ax_l.set_title('longitudinal C(z)')
    ax_l.set_ylabel(r"$\Delta\phi$ [deg]")
    ax_l.tick_params(axis='y', colors='green')
    ax_l.yaxis.label.set_color('green')
    ax_l.plot(z,cz,label=r"$\Delta\phi$",color='green')
    # ax_l.plot(z,sz,color='green')
    vscale=ax_l.axis()[3]*0.1
    viseoz = [x*vscale for x in viseo]
    ax_l.plot(z,viseoz,label='',color='black')
    ax_l.plot(z,zero,color='black')

    ax_r = ax_l.twinx()
    ax_r.set_ylabel(r'$\Delta$w/w [%]')
    ax_r.tick_params(axis='y', colors='red')
    ax_r.yaxis.label.set_color('red')
    ax_r.plot(z,cdw,label=r'$\Delta$w/w',color='red')
    # ax_r.plot(z,sdw,color='red')
    ax_r.plot(z,zero,color='red', linestyle='--')
    #----------*----------*
    show()

def loesung(filepath):                 ## total classic FODO lattice
    lattice = parse_yaml_and_fabric(filepath)
    soll_track = track_soll(lattice)   ## !WICHTIG! track soll Teilchen hier
    lattice.stats(soll_track)          ## count elements and other statistics
    epsiz()                            ## longitudinal params
    ## Rechne: ganze Zelle und Anfangswerte
    mcell,betax,betay = lattice.cell(closed=CONF['periodic'])
    collect_summary()
    dictprnt(SUMMARY,text='summary')
    ## Zeige Grafik: Lösungen als Funktion von (s)
    functions = lattice.twiss_functions(30)
    display(functions)   # twiss functions
## main-------
if __name__ == '__main__':
    import sys
    # filepath = 'LEBT_fodo_with_RF.yml'       ## input file (syntax=YAML)
    filepath = 'LEBT_fodo_with_RF(x).yml'      ## input file (syntax=YAML)
    if len(sys.argv) == 2:
        filepath = sys.argv[1]
    loesung(filepath)
