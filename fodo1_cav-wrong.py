# -*- coding: utf-8 -*-
from setup import Phys
from elements import k0,I,D,QF,QD,SD,WD,CAV,RFG,Beam
from lattice import Lattice
from pylab import plot, show, legend
from math import sqrt

#Werte={'lqd','lqf','ld','lcav','U0','phi0','fRF','dBdz','dWf','verbose'}
Werte ={} # Eigabewerte fuer eine basis zelle (als gobal definiert! ..bad but efficient.)

def plotter(beta_fun,cos_like,sin_like):  ## plotting
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
def make_cavity(l):   ## kavität
    global Werte
    w = Werte
    tk = Beam.soll.tkin                    # kinetic energy @ entrance
    cavity = Lattice()
    dri = D(length=0.5*l,beam=Beam.soll)   # drift before RFgap
    gap = CAV(U0=w['U0'],PhiSoll=w['phi0'],fRF=w['fRF'],label='cav',beam=Beam.soll,dWf=w['dWf'])  # T.Wrangler, Dr.Tiede
    # gap = RFG(U0=w['U0'],PhiSoll=w['phi0'],fRF=w['fRF'],label='rfg',beam=Beam.soll,dWf=w['dWf'])  # Trace3D
    drf = D(length=0.5*l,beam=Beam.soll)   # drift after RFgap
    cavity.add_element(dri)
    cavity.add_element(gap)
    cavity.add_element(drf)
    # dictp(Beam.soll,'gap exit',{'matrix'})
    return cavity
def make_rf_section(lcav,gaps=1):   ## RF sektion
    ''' gaps: nboff gaps per rf section'''
    global Werte
    l = 0.                          # length counter
    section = Lattice()
    for i in range(gaps):
        cav = make_cavity(lcav)
        section.append(cav)
        l += lcav
    return section  
def make_half_cell(upstream=True,gaps=3):  # 1/2 cell
    global Werte
    w = Werte
    tki = Beam.soll.tkin                            # kinetic energy of particle @ entrance
    # basis zelle
    cell=Lattice()
    if upstream : # 1/2 basis zelle upstream
        tk = tki                                     
        kq   = k0(gradient=w['dBdz'],tkin=tk)                        # quad strength @ entrance
        mqf=QF(k0=kq,length=w['lqf'],label='QF',beam=Beam.soll)      # F quad before cavities
        rf_section = make_rf_section(w['lcav'],gaps)                 # cavities
        ld = w['ld']-rf_section.length                               # fit drift glue
        md=D(0.5*ld,beam=Beam.soll)                                  # drift zw. F quad und cavities
        cell.add_element(mqf)
        cell.add_element(md)
        cell.append(rf_section)

        tk = Beam.soll.tkin                                          # new energy
        kq   = k0(gradient=w['dBdz'],tkin=tk)                        # new quad strength 
        mqd=QD(k0=kq,length=w['lqd'],label='QD',beam=Beam.soll)      # D quad after cavities
        md=D(0.5*ld,beam=Beam.soll)                                  # drift zw. cavities und D quad
        cell.add_element(md)
        cell.add_element(mqd)
        # cell.out()
    else:  #  1/2 basis zelle downstream  = reverse of upstream
        tk = tki               
        kq   = k0(gradient=w['dBdz'],tkin=tk) 
        mqd=QD(k0=kq,length=w['lqd'],label='QD',beam=Beam.soll)      
        rf_section = make_rf_section(w['lcav'],gaps)  
        ld = w['ld']-rf_section.length
        md=D(0.5*ld,beam=Beam.soll)                   
        cell.add_element(mqd)
        cell.add_element(md)
        cell.append(rf_section)
 
        tk = Beam.soll.tkin
        kq   = k0(gradient=w['dBdz'],tkin=tk)    
        mqf=QF(k0=kq,length=w['lqf'],label='QF',beam=Beam.soll) 
        md=D(0.5*ld,beam=Beam.soll)                  
        cell.add_element(md)
        cell.add_element(mqf)
        # cell.out()
    deltaTK=Beam.soll.tkin - tki
    return cell,deltaTK
#############################################################################
def dictp(what,text='========',filter={}):  ## helper to print dicts
        print('========= '+text+' =================')
        for k,v in what.__dict__.items():
            if k in filter:
                continue
            print(k.rjust(30),':',v)
def loesung1():  # total classic FODO lattice
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
    
    sollParticle = Beam.soll                 # sollParticlle @ injection energy
    gaps_per_half_cell= 3                    # KNOB:  gaps/cell

    dWf=1.0                                  # acceleration                   
    dBdz0  = Phys['quad_gradient']*7.85      # KNOB: quad gradient
    dBdz0  = Phys['quad_gradient']*9.9375      # KNOB: quad gradient
    nboff_super_cells = 15*10                # KNOB:  final energy
    nboff_super_cells = 15*8                 # KNOB:  final energy
    # nboff_super_cells = 15*5                 # KNOB:  final energy
    # nboff_super_cells = 15*1                 # KNOB:  final energy

    # dWf=0.                                   # no acceleration
    # dBdz0  = Phys['quad_gradient']*9.        # KNOB: quad gradient
    # nboff_super_cells = 15*1                 # KNOB:  final energy
    # nboff_super_cells = 15*8                 # KNOB:  final energy

    global Werte           # store globals
    Werte={'lqd':lqd,'lqf':lqf,'ld':ld,'lcav':lcav,'U0':u0,'phi0':phi0,'fRF':fRF0,'dBdz':dBdz0,'dWf':dWf,'verbose':True}
    #-----------------------------------------
    super_cell=Lattice()
    nboff_gaps=0                 # gap counter
    tk = tk0                     # energy counter
    for icell in range(nboff_super_cells):
        #------
        # kann man die struktur bei höheren energien ändern?
        # if tk >= 150.:                # KNOB: energy at which...
            # gaps_per_half_cell=3                 # KNOB: change gaps/cell
        #------
        zelle = Lattice()  # basis zelle
        (half_cell,deltaTK) = make_half_cell(upstream=True,gaps=gaps_per_half_cell); nboff_gaps+=gaps_per_half_cell
        zelle.append(half_cell)
        tk += deltaTK
        (half_cell,deltaTK) = make_half_cell(upstream=False,gaps=gaps_per_half_cell); nboff_gaps+=gaps_per_half_cell
        zelle.append(half_cell)
        tk += deltaTK  # new energy update here!
        # cell.out()
        super_cell.append(zelle)  # add zelle to super cell
    # super_cell.out()
    lattice_length=super_cell.length
    # print('lattice length [m]={}'.format(lattice_length))
    #-----------------------------------------
    # Anfangswerte
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