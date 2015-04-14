# -*- coding: utf-8 -*-
import setup as IN
from elements import I,D,QF,QD,SD,WD,CAV
from lattice import Lattice
from pylab import plot, show, legend
from math import sqrt

def plotter(beta_fun,cos_like,sin_like):
    s   = [ x[0] for x in beta_fun]    # bahnkoordinate s
    bx  = [ x[1] for x in beta_fun]    # beta-x
    bxn = [-x[1] for x in beta_fun]    # beta-x (negatif)
    by  = [ x[2] for x in beta_fun]    # beta-y
    byn = [-x[2] for x in beta_fun]    # beta-y (negatif)
    
    cx = [x[0] for x in cos_like]   # cos-like-x
    cy = [x[2] for x in cos_like]   # cos-like-y
    sx = [x[0] for x in sin_like]   # sin-like-x
    sy = [x[2] for x in sin_like]   # sin-like-x
    
    viseo = [x[3] for x in beta_fun]
    zero=[0. for x in beta_fun]
    
    plot(s,bx ,label='betax',color='green')
    # plot(s,bxn,label='',     color='green')
    plot(s,by ,label='betay',color='red')
    # plot(s,byn,label='',     color='red')
    
    plot(s,cx,label='Cx(s)',color='blue')
    plot(s,sx,label='Sx(s)',color='brown') 
    # plot(s,cy,label='Cy(s)',color='blue')
    # plot(s,sy,label='Sy(s)',color='brown')
    
    plot(s,viseo,label='',color='black')
    plot(s,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    show()
def loesung1():
    # physik werte
    tk    = IN.physics['kinetic_energy']
    Bgrad = IN.physics['quad_gradient']
    k0    = IN.k0_p(gradient=Bgrad,tkin=tk)  # quad strength
    print('k0[1/m**2]= ',k0)
    kqf= k0*8.8
    kqd= kqf 
    
    # längen
    lqd=  0.4       # QD len
    lqf=  0.5*lqd   # QF len
    # focal = 1./(kqd*lqd)
    # ld=2.*sqrt(2.)*focal    # from Roßbach CERN-94-01-V1 pp.65
    ld   =  0.4     # drift len
    lcav = 0.08     # cav len
    print('l-quad = ',lqd)
    print('l-drift= ',ld)
    print('l-cav  = ',lcav)
    print()

    # elemente
    mqf=QF(kqf,lqf,'QF') 
    mqd=QD(kqd,lqd,'QD')
    mcl = mcr = D(length=0.5*lcav,label='cav')
    cavity=CAV(
    # __init__(self, U0=10., TrTF=0.5, PhiSoll=-0.25*pi, Tkin=50., fRF=800., label='CAV'):
        # U0=IN.physics['spalt_spannung'],
        U0=3.0,
        PhiSoll=IN.physics['soll_phase']*IN.physics['radians'],
        Tkin=IN.physics['kinetic_energy'],
        fRF=IN.physics['frequenz'],
        label='gap')
    
    # kavität
    cav=Lattice()
    cav.add_element(mcr)
    cav.add_element(cavity)
    cav.add_element(mcl)
    # cav.out()
    
    # RF sektion
    cnt1=0   # cav/section
    rf_section = Lattice()
    rf_section.append(cav); cnt1+=1
    rf_section.append(cav); cnt1+=1
    rf_section.append(cav); cnt1+=1
    # rf_section.out()  
    
    # abgleich drift strecke
    lrf_section=rf_section.length
    ld=ld-lrf_section
    print('l-drift = ',ld)
    print('l-cav_section  = ',lrf_section)
    print()
    md=D(0.5*ld)

    # basis zelle
    cnt2=0   # sec/cell
    cell=Lattice()
    cell.add_element(mqf)
    cell.add_element(md)
    cell.append(rf_section);cnt2+=1
    cell.add_element(md)
    cell.add_element(mqd)
    cell.add_element(md)
    cell.append(rf_section);cnt2+=1
    cell.add_element(md)
    cell.add_element(mqf)
    cell_length=cell.length
    # cell.out()
    
    # mehrere zellen (several cells)
    cnt3=0   # cells/super_cell
    super_cell=Lattice()
    super_cell.append(cell);cnt3+=1
    super_cell.append(cell);cnt3+=1
    super_cell.append(cell);cnt3+=1
    super_cell.append(cell);cnt3+=1
    lattice_length=super_cell.length
    print('lattice length [m]={}'.format(lattice_length))
    # super_cell.out()
    
    # anfangswerte
    mcell,betax,betay=super_cell.cell()
    print()
    print('BETAx[0] {:.3f} BETAy[0] {:.3f}'.format(betax,betay))
    
    # lösungen als f(s)
    beta_func = super_cell.beta_functions(20)   
    cossin_like = super_cell.cossin(20)
    
    # Zusammenfassung
    s_p=IN.Proton(tk)
    s_name = s_p.name
    s_e0 = s_p.e0
    s_gaps=cnt1*cnt2*cnt3
    s_bgrad=IN.dBdz_p(kqf,tk)
    s_u0=cavity.u0
    s_utot=s_u0*s_gaps
    s_accel=s_utot/lattice_length
    s_l200= 200./s_accel
    summary={
    'quadrupole size          [m]':lqd,
    'particle rest mass[MeV/c**2]':s_e0,
    'particle energy        [Mev]':tk,
    'qudrupole strength  [1/m**2]':kqf,
    'number of cavities          ':s_gaps,
    'qudrupole gradient     [T/m]':s_bgrad,
    'av. voltage_per_gap     [MV]':s_u0,
    'lattice_length           [m]':lattice_length,
    'tot.acceleration       [MeV]':s_utot,
    'av.acceleration      [MeV/m]':s_accel,
    'accl. length for 200 Mev [m]':s_l200,
    }
    print('\nSummary('+s_name+')')
    for s in summary.items():
        print('{}=  {:.3f}'.format(s[0].rjust(44),s[1]))
    
    # grafik
    plotter(beta_func,cossin_like[0],cossin_like[1])
if __name__ == '__main__':
    loesung1()