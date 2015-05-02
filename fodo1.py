#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import Phys
from elements import k0,I,D,QF,QD,SD,WD,CAV,RFG,Beam
from lattice import Lattice
from pylab import plot, show, legend, figure, subplot, axis
from math import sqrt

#Werte={'lqd','lqf','ld','lcav','U0','phi0','fRF','dBdz','dWf','verbose'}
Werte ={} # Eigabewerte fuer eine basis zelle (als gobal definiert! ..bad but efficient.)

def display(functions):  ## plotting
    #----------*----------*   # unpack
    beta_fun = functions[0]
    cos_like = functions[1]
    sin_like = functions[2]
    emix=Phys['emitx(i)']  # emittance @ entrance
    emiy=Phys['emity(i)']  # emittance @ entrance
    #----------*----------*   # bahnkoordinate z
    z   = [ x[0] for x in beta_fun]    
    #----------*----------*
    bx  = [ sqrt(x[1]*emix) for x in beta_fun]    # envelope (beta-x)
    by  = [ sqrt(x[2]*emiy) for x in beta_fun]    # envelope (beta-y)
    bxn = [-x for x in bx]    # beta-x (negatif)
    byn = [-x for x in by]    # beta-y (negatif)
    zero= [0. for x in beta_fun]  # zero line
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
    width=20; height=12
    figure('FODO 1',figsize=(width,height))
    #----------*----------*   # transverse X
    splot=subplot(311)
    splot.set_title('transverse x')
    plot(z,bx ,label=r'$\sigma$ [m]',color='green')
    plot(z,bxn,label='',color='green')
    plot(z,cx,label='Cx[m]',color='blue',linestyle='-.')
    plot(z,sx,label='Sx[m]',color='red' ,linestyle='-.') 
    vscale=axis()[3]*0.1
    viseo = [x[3]*vscale for x in beta_fun]
    plot(z,viseo,label='',color='black')
    plot(z,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    #----------*----------*   # transverse Y
    splot=subplot(312)
    splot.set_title('transverse y')
    plot(z,by ,label=r'$\sigma$ [m]',color='green')
    plot(z,byn,label='',color='green')
    plot(z,cy,label='Cy[m]',color='blue',linestyle='-.')
    plot(z,sy,label='Sy[m]',color='red' ,linestyle='-.')
    vscale=axis()[3]*0.1
    viseo = [x[3]*vscale for x in beta_fun]
    plot(z,viseo,label='',color='black')
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
    viseo = [x[3]*vscale for x in beta_fun]
    ax_l.plot(z,viseo,label='',color='black')
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
def make_cavity(l):   ## kavität
    global Werte
    w = Werte
    tk = Beam.soll.tkin                    # kinetic energy @ entrance
    cavity = Lattice()
    dri = D(length=0.5*l,beam=Beam.soll,label='>')   # drift before RFgap
    # gap = CAV(U0=w['U0'],PhiSoll=w['phi0'],fRF=w['fRF'],label='cav',beam=Beam.soll,dWf=w['dWf'])  # T.Wrangler, Dr.Tiede
    gap = RFG(U0=w['U0'],PhiSoll=w['phi0'],fRF=w['fRF'],label='rfg',beam=Beam.soll,dWf=w['dWf'])  # Trace3D
    drf = D(length=0.5*l,beam=Beam.soll,label='<')   # drift after RFgap
    cavity.add_element(dri)
    cavity.add_element(gap)
    cavity.add_element(drf)
    Phys['LCAV']=cavity.length
    # cavity.out()
    # objprnt(Beam.soll,'gap exit',{'matrix'})
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
    Phys['RFSection']=section.length
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
    Phys['LQF']=2.*mqf.length
    Phys['LQD']=2.*mqd.length
    Phys['LD'] =md.length
    Phys['CELL']=cell.length
    deltaTK=Beam.soll.tkin - tki
    return cell,deltaTK
#############################################################################
def objprnt(what,text='========',filter={}):  ## helper to print objects as dictionary
        print('\n========= '+text+' =================')
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
    
    gaps_per_half_cell= 3                    # KNOB:  gaps/cell

    dBdz0  = Phys['quad_gradient']*7.85      # KNOB: quad gradient
    # nboff_super_cells = 15*10                # KNOB:  final energy
    nboff_super_cells = 15*8                 # KNOB:  final energy
    # nboff_super_cells = 15*5                 # KNOB:  final energy
    nboff_super_cells = 15*1                 # KNOB:  final energy
    # nboff_super_cells = 1                    # KNOB:  final energy

    dWf=1.0                                  # acceleration flag=yes              
    # dWf=0.                                   # acceleration flag=no
    # dBdz0  = Phys['quad_gradient']*9.        # KNOB: quad gradient
    # nboff_super_cells = 15*1                 # KNOB:  final energy
    # nboff_super_cells = 15*8                 # KNOB:  final energy

    global Werte           # store globals
    Werte={'lqd':lqd,'lqf':lqf,'ld':ld,'lcav':lcav,'U0':u0,'phi0':phi0,'fRF':fRF0,'dBdz':dBdz0,'dWf':dWf,'verbose':True}
    w = Werte
    #--------------------1st cavity-----------
    gapi = RFG(U0=w['U0'],PhiSoll=w['phi0'],fRF=w['fRF'],label='rfg',beam=Beam.soll,dWf=w['dWf'])  # Trace3D
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
        if icell == 0:
            base_cell = zelle
        super_cell.append(zelle)  # add zelle to super cell
    # super_cell.out()
    lattice_length=super_cell.length
    # print('lattice length [m]={}'.format(lattice_length))
    #-----------------------------------------
    # Berechne ganze Zelle und Anfangswerte 
    mcell,betax,betay=super_cell.cell(closed=False)
    print('\nBETAx(i) {:.3g} [m], BETAy(i) {:.3g} [m]'.format(betax,betay))
    #---------------last cavity---------------
    gapf = RFG(U0=w['U0'],PhiSoll=w['phi0'],fRF=w['fRF'],label='rfg',beam=Beam.soll,dWf=w['dWf'])  # Trace3D
    # objprnt(gapi, filter={'matrix'})
    # objprnt(gapf, filter={'matrix'})    
    #-----------------------------------------
    # Zusammenfassung
    s_ttf_i =gapi.tr
    s_ttf_f =gapf.tr
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
    s_emi   =Phys['emitx(i)']
    s_aper  =Phys['sigx(i)']
    s_phis  =Phys['soll_phase']
    s_lamb  =Phys['wellenlänge']
    s_freq  =Phys['frequenz']
    s_nboff_cells =nboff_super_cells
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
    'emittance (i)        [m*rad]':s_emi,
    'sigma-x,-y (i)           [m]':s_aper,
    'sync. phase            [deg]':s_phis,
    'wave length              [m]':s_lamb,
    'frequency              [MHz]':s_freq,
    'time transition factor (i,f)':(s_ttf_i,s_ttf_f),
    'number of cells             ': s_nboff_cells,
    }
    print('\n============= Summary =================')
    for k,v in summary.items():
        print(k.rjust(30),':',v)
    print('\n====== base cell ========:\n'+zelle.report())
    #-----------------------------------------
    # Grafik: lösungen als f(s)
    functions = super_cell.functions(30)   
    display(functions)
if __name__ == '__main__':
    loesung1()