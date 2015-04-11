# -*- coding: utf-8 -*-
import setup as IN
from elements import I,D,QF,QD,SD,WD,CAV
from lattice import Lattice
from pylab import plot, show, legend

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
    plot(s,bxn,label='',     color='green')
    plot(s,by ,label='betay',color='red')
    plot(s,byn,label='',     color='red')
    
    plot(s,cx,label='Cx(s)',color='blue')
    plot(s,sx,label='Sx(s)',color='brown') 
    plot(s,cy,label='Cy(s)',color='blue')
    plot(s,sy,label='Sy(s)',color='brown')
    
    plot(s,viseo,label='',color='black')
    plot(s,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    show()

def test1():
    # physik werte
    tk = IN.physics['kinetic_energy']
    Bgrad = IN.physics['quad_gradient']
    k0 = IN.k0p(gradient=Bgrad,tkin=tk)  # quad strength
    print('k0[1/m**2]= ',k0)
    kqf= k0*2.5
    kqd= kqf
    
    # längen
    lqd=  0.4       # QD len
    lqf=  0.5*lqd   # QF len
    ld =  0.8       # drift len
    lcav = 0.08     # cav len

    # elemente
    mqf=QF(kqf,lqf,'QF') 
    mqd=QD(kqd,lqd,'QD')
    mcl = mcr = D(length=0.5*lcav,label='cav')
    # __init__(self, U0=10., TrTF=0.5, PhiSoll=-0.25*pi, Tkin=50., fRF=800., label='CAV'):
    cavity=CAV(
        # U0=IN.physics['spalt_spannung'],
        U0=1.,
        TrTF=IN.physics['transit_time'],
        PhiSoll=IN.physics['soll_phase']*IN.physics['radians'],
        Tkin=IN.physics['kinetic_energy'],
        fRF=IN.physics['frequenz'],
        label='gap')
    
    # the cavity
    cav=Lattice()
    cav.add_element(mcl)
    cav.add_element(cavity)
    cav.add_element(mcr)
    
    #rf section
    count = 0
    rf_section = Lattice()
    rf_section.append(cav); count +=1
    rf_section.append(cav); count +=1
    rf_section.append(cav); count +=1
    rf_section.append(cav); count +=1
    rf_section.append(cav); count +=1
    rf_section.append(cav); count +=1
    rf_section.append(cav); count +=1
    rf_section.append(cav); count +=1
    
    # drift glue
    md=D(0.5*(ld-count*lcav))

    ## lattice
    lattice=Lattice()
    lattice.add_element(mqf)
    lattice.add_element(md)
    lattice.append(rf_section)
    lattice.add_element(md)
    lattice.add_element(mqd)
    lattice.add_element(md)
    lattice.append(rf_section)
    lattice.add_element(md)
    lattice.add_element(mqf)
    
    lattice.append(lattice)
    # lattice.append(lattice)
    # lattice.append(lattice)
    # lattice.append(lattice)
    # lattice.out()
    ## cell boundaries
    mcell,betax,betay=lattice.cell()
    print('BETAx[0] {:.3f} BETAy[0] {:.3f}'.format(betax,betay))
    ## lattice function as f(s)
    beta_func = lattice.beta_functions(120)   
    cossin_like = lattice.cossin(120)
    ## plots
    plotter(beta_func,cossin_like[0],cossin_like[1])
####################################################
if __name__ == '__main__':
    test1()