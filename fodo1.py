# -*- coding: utf-8 -*-
from setup import Beam, Phys
from elements import k0,I,D,QF,QD,SD,WD,CAV,RFG
from lattice import Lattice
from pylab import plot, show, legend
from math import sqrt

#Werte={'lqd','lqf','ld','lcav','U0','phi0','fRF','dBdz','dWf','verbose'}
Werte ={} # Eigabewerte fuer eine basis zelle (als gobal definiert! ..bad but efficient.)

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
    plot(s,by ,label='betay',color='red')
    # plot(s,bxn,label='',     color='green')
    # plot(s,byn,label='',     color='red')
    
    plot(s,cx,label='Cx(s)',color='blue')
    plot(s,sx,label='Sx(s)',color='brown') 
    # plot(s,cy,label='Cy(s)',color='green')
    # plot(s,sy,label='Sy(s)',color='red')
    
    plot(s,viseo,label='',color='black')
    plot(s,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    show()
def make_cavity(tki,l):   # kavität
    global Werte
    w = Werte
    tk = tki        # kinetic energy @ entrance
    beam = Beam(tk)
    cavity = Lattice()
    dri = D(length=0.5*l,beam=beam)   # drift before RFgap
    gap = RFG(U0=w['U0'],PhiSoll=w['phi0'],fRF=w['fRF'],label='rfg',beam=beam,dWf=w['dWf'])
    # if w['verbose']':
        # print('========= GAP =================')
        # for k,v in gap.__dict__.items():
            # if k=='matrix':
                # continue
            # print(k.rjust(30),':',v)
    tk += gap.deltaW        # kinetic energy after RFgap
    beam = Beam(tk)
    drf = D(length=0.5*l,beam=beam)     # drift after RFgap
    cavity.add_element(dri)
    cavity.add_element(gap)
    cavity.add_element(drf)
    return (cavity,tk)
def make_rf_section(tki,lcav,gaps=1):   # RF sektion
    ''' gaps: nboff gaps per rf section'''
    tk = tki        # kinetic energy @ entrance
    l = 0.          # length of section
    section = Lattice()
    for i in range(gaps):
        (cav,tkf) = make_cavity(tk,lcav)
        section.append(cav)
        l += lcav
        tk = tkf    # kinetic energy after cavity
    return (section,tkf)  
def make_half_cell(tki,upstream=True,gaps=3):
    global Werte
    w = Werte
    # basis zelle
    cell=Lattice()
    if upstream : # 1/2 basis zelle upstream
        tk   = tki                              # kinetic energy @ entrance
        kq   = k0(gradient=w['dBdz'],tkin=tk)   # quad strength @ entrance
        beam = Beam(tk)                         # beam @ entrance
        mqf=QF(k0=kq,length=w['lqf'],label='QF',beam=beam)      # F quad before cavities
        (rf_section,tkf) = make_rf_section(tk,w['lcav'],gaps)   # cavities
        ld = w['ld']-rf_section.length
        md=D(0.5*ld,beam=beam)                   # drift zw. F quad und cavities
        cell.add_element(mqf)
        cell.add_element(md)
        cell.append(rf_section)
        tk   = tkf                               # energy update after cavities
        kq   = k0(gradient=w['dBdz'],tkin=tk)    # new quad strength 
        beam = Beam(tk)                          # new beam
        mqd=QD(k0=kq,length=w['lqd'],label='QD',beam=beam) # D quad after cavities
        md=D(0.5*ld,beam=beam)                   # drift zw. cavities und D quad
        cell.add_element(md)
        cell.add_element(mqd)
        # cell.out()
    else:  #  1/2 basis zelle downstream  = reverse of upstream
        tk   = tki
        kq   = k0(gradient=w['dBdz'],tkin=tk)        
        beam = Beam(tk)                             
        mqd=QD(k0=kq,length=w['lqf'],label='QF',beam=beam) 
        (rf_section,tkf) = make_rf_section(tk,w['lcav'],gaps)
        ld = w['ld']-rf_section.length
        md=D(0.5*ld,beam=beam)
        cell.add_element(mqd)
        cell.add_element(md)
        cell.append(rf_section)
        tk   = tkf
        kq   = k0(gradient=w['dBdz'],tkin=tk)           
        beam = Beam(tk)
        mqf=QF(k0=kq,length=w['lqd'],label='QD',beam=beam) 
        md=D(0.5*ld,beam=beam) 
        cell.add_element(md)
        cell.add_element(mqf)
        # cell.out()
    deltaTK=tkf-tki
    return cell,deltaTK
#############################################################################
def loesung1():
    #-----------------------------------------
    # längen
    lqd  =  0.2     # 1/2 QD len
    lqf  =  0.2     # 1/2 QF len
    ld   =  0.4     # drift len            # KNOB: effective focus of FODO
    lcav =  0.08    # cav len
    #-----------------------------------------
    # physik werte
    u0     = Phys['spalt_spannung']
    phi0   = Phys['soll_phase']*Phys['radians']
    fRF0   = Phys['frequenz']
    tk0    = Phys['kinetic_energy']*1.       # KNOB: injection energy
    dBdz0  = Phys['quad_gradient']*8.05      # KNOB: quad gradient
    dBdz0  = Phys['quad_gradient']*7.85      # KNOB: quad gradient
    global Werte            # store globals
    Werte={'lqd':lqd,'lqf':lqf,'ld':ld,'lcav':lcav,'U0':u0,'phi0':phi0,'fRF':fRF0,'dBdz':dBdz0,'dWf':1.,'verbose':False}
    #-----------------------------------------
    super_cell=Lattice()
    nboff_gaps=0                 # gap counter
    tk = tk0                     # energy counter
    # nboff_super_cells = 15*5   # KNOB:  final energy
    nboff_super_cells = 15*8     # KNOB:  final energy
    # nboff_super_cells = 15*1   # KNOB:  final energy
    gaps_per_half_cell= 3        # KNOB:  gaps/cell
    for icell in range(nboff_super_cells):
        #------
        # kann man die struktur bei höheren energien ändern?
        # if tk >= 150.:                # KNOB: energy at which...
            # gaps_per_half_cell=3                 # KNOB: change gaps/cell
        #------
        zelle = Lattice()  # basis zelle
        (half_cell,deltaTK) = make_half_cell(tk,upstream=True,gaps=gaps_per_half_cell); nboff_gaps+=gaps_per_half_cell
        zelle.append(half_cell)
        tk += deltaTK
        (half_cell,deltaTK) = make_half_cell(tk,upstream=False,gaps=gaps_per_half_cell); nboff_gaps+=gaps_per_half_cell
        zelle.append(half_cell)
        tk += deltaTK  # new energy update here!
        # cell.out()
        super_cell.append(zelle)  # add zelle to super cell
    # super_cell.out()
    lattice_length=super_cell.length
    # print('lattice length [m]={}'.format(lattice_length))
    #-----------------------------------------
    # anfangswerte
    mcell,betax,betay=super_cell.cell(closed=True)
    print()
    print('BETAx[0] {:.3f} BETAy[0] {:.3f}'.format(betax,betay))
    #-----------------------------------------
    # Zusammenfassung
    s_tk_i  =tk0
    s_tk_f  =tk
    s_lqd   =lqd
    s_p     =Beam(s_tk_f)
    s_name  =s_p.name
    s_e0    =s_p.e0
    s_gaps  =nboff_gaps
    s_bgrad =dBdz0
    s_kq_i  =k0(gradient=s_bgrad,tkin=s_tk_i)
    s_kq_f  =k0(gradient=s_bgrad,tkin=s_tk_f)
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