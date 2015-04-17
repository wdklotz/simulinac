# -*- coding: utf-8 -*-
import setup as UTIL
from elements import I,D,QF,QD,SD,WD,CAV
from lattice import Lattice
from pylab import plot, show, legend
from math import sqrt
#Werte={'lqd':lqd,'lqf':lqf,'ld':ld,'lcav':lcav,'U0':u0,'phi0':phi0,'fRF':fRF0,'tkin':tk0,'dBdz':dBdz0}
Werte ={} # Eigabewerte fuer eine basis zelle
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
    
    # plot(s,cx,label='Cx(s)',color='blue')
    plot(s,sx,label='Sx(s)',color='brown') 
    plot(s,cy,label='Cy(s)',color='blue')
    # plot(s,sy,label='Sy(s)',color='brown')
    
    plot(s,viseo,label='',color='black')
    plot(s,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    show()
def make_half_cell(upstream=True,verbose=False):
    global Werte
    w = Werte
    #-----------------------------------------
    # elemente
    tkin = w['tkin']                                # updated kinetic energy
    kq=UTIL.k0(gradient=w['dBdz'],tkin=tkin)        # update quad strength
    mqf=QF(kq,w['lqf'],'QF')                        # update F quad
    mqd=QD(kq,w['lqd'],'QD')                        # update D quad
    cavity = CAV(U0=w['U0'],PhiSoll=w['phi0'],Tkin=tkin,fRF=w['fRF'],label='gap')  # update cavity
    mcl = mcr = D(length=0.5*w['lcav'],label='cav') # drifts do not change
    if verbose:
        print('========= CAVITY =================')
        for k,v in cavity.__dict__.items():
            if k=='matrix' or k=='prot':
                continue
            print(k.rjust(30),':',v)
        print('========= QUADRUPOLE ==============')
        for k,v in mqd.__dict__.items():
            if k=='matrix':
                continue
            print(k.rjust(30),':',v)
    #-----------------------------------------
    # kavität
    cav=Lattice()
    cav.add_element(mcr)
    cav.add_element(cavity)
    cav.add_element(mcl)
    # cav.out()
    #-----------------------------------------
    # RF sektion
    cnt1=0   # cav/section
    rf_section = Lattice()
    rf_section.append(cav); cnt1+=1
    rf_section.append(cav); cnt1+=1
    rf_section.append(cav); cnt1+=1
    # rf_section.out()  
    #-----------------------------------------
    # abgleich drift strecke
    lrf_section=rf_section.length
    ld=w['ld']-lrf_section
    # print('Abgleich')
    # print('l-drift = ',ld)
    # print('l-cav_section  = ',lrf_section)
    # print()
    md=D(0.5*ld) # drift zw. cavity u. quad
    #-----------------------------------------
    # basis zelle
    cell=Lattice()
    cnt2=0
    if upstream : # 1/2 basis zelle upstream
        cell.add_element(mqf)
        cell.add_element(md)
        cell.append(rf_section);cnt2+=1
        cell.add_element(md)
        cell.add_element(mqd)
        # cell.out()
    else:  #  1/2 basis zelle downstream  
        cell.add_element(mqd)
        cell.add_element(md)
        cell.append(rf_section);cnt2+=1
        cell.add_element(md)
        cell.add_element(mqf)
    nboff_gaps = cnt1*cnt2
    dW=nboff_gaps*cavity.deltaW
    return cell,nboff_gaps,dW
def loesung1():
    global Werte
    #-----------------------------------------
    # längen
    lqd  =  0.2     # 1/2 QD len
    lqf  =  0.2     # 1/2 QF len
    ld   =  0.4     # drift len
    lcav =  0.08    # cav len
    # physik werte
    #-----------------------------------------
    u0     = UTIL.physics['spalt_spannung']
    phi0   = UTIL.physics['soll_phase']*UTIL.physics['radians']
    fRF0   = UTIL.physics['frequenz']
    tk0    = UTIL.physics['kinetic_energy']
    dBdz0  = UTIL.physics['quad_gradient']*8.5
    #-----------------------------------------
    super_cell=Lattice()
    nboff_super_cells = 16
    nboff_gaps=0
    Werte={'lqd':lqd,'lqf':lqf,'ld':ld,'lcav':lcav,'U0':u0,'phi0':phi0,'fRF':fRF0,'tkin':tk0,'dBdz':dBdz0}
    for icell in range(nboff_super_cells):
        # basis zelle
        cell = Lattice()   # updated cell
        (half_cell,cnt,deltaW)   = make_half_cell(upstream=True); nboff_gaps+=cnt
        cell.append(half_cell)
        Werte['tkin'] += deltaW  # energy update here!
        (half_cell,cnt,deltaW) = make_half_cell(upstream=False); nboff_gaps+=cnt
        cell.append(half_cell)
        Werte['tkin'] += deltaW  # energy update here!
        # cell.out()
        super_cell.append(cell)  # add to super cell
    # super_cell.out()
    lattice_length=super_cell.length
    # print('lattice length [m]={}'.format(lattice_length))
    #-----------------------------------------
    # anfangswerte
    mcell,betax,betay=super_cell.cell()
    print()
    print('BETAx[0] {:.3f} BETAy[0] {:.3f}'.format(betax,betay))
    #-----------------------------------------
    # Zusammenfassung
    s_tk_i  =tk0
    s_tk_f  =Werte['tkin']
    s_lqd   =lqd
    s_p     =UTIL.Proton(s_tk_f)
    s_name  =s_p.name
    s_e0    =s_p.e0
    s_gaps  =nboff_gaps
    s_bgrad =dBdz0
    s_kq_i  =UTIL.k0(gradient=s_bgrad,tkin=s_tk_i)
    s_kq_f  =UTIL.k0(gradient=s_bgrad,tkin=s_tk_f)
    s_u0    =u0
    s_utot  =s_tk_f - s_tk_i
    s_latlen=lattice_length
    s_accel =s_utot/s_latlen
    summary={
    'quadrupole size          [m]':s_lqd,
    'particle rest mass[MeV/c**2]':s_e0,
    'particle energy(i)     [Mev]':s_tk_i,
    'particle energy(f)     [Mev]':s_tk_f,
    'quad strength(i)    [1/m**2]':s_kq_i,
    'quad strength(f)    [1/m**2]':s_kq_f,
    'number of cavities          ':s_gaps,
    'qudrupole gradient     [T/m]':s_bgrad,
    'av. voltage_per_gap     [MV]':s_u0,
    'lattice_length           [m]':s_latlen,
    'tot.acceleration       [MeV]':s_utot,
    'av.acceleration      [MeV/m]':s_accel,
    }
    print('\n======== Summary('+s_name+') =========')
    for k,v in summary.items():
        print(k.rjust(30),':',v)
    #-----------------------------------------
    # Grafik: lösungen als f(s)
    beta_func   = super_cell.beta_functions(20)   
    cossin_like = super_cell.cossin(20)
    plotter(beta_func,cossin_like[0],cossin_like[1])
if __name__ == '__main__':
    loesung1()