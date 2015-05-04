#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import Phys,k0,dictprnt,Beam
from elements import D,QF,QD,RFG,CAV,RFC
from lattice import Lattice
from pylab import plot,show,legend,figure,subplot,axis
from math import sqrt

def display(functions):  ## plotting
    #----------*----------*   # unpack
    beta_fun = functions[0]
    cos_like = functions[1]
    sin_like = functions[2]
    emix=Phys['emitx_i']  # emittance @ entrance
    emiy=Phys['emity_i']  # emittance @ entrance
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
    # figure('FODO 1')
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
def make_rf_section(w):   ## RF sektion
    gaps = w['gaps']    # gaps/half-cell
    section = Lattice()
    for i in range(gaps):
        cav=RFC(length=w['lcav'],U0=w['U0'],PhiSoll=w['phi0'],fRF=w['fRF'],label='rfg',beam=Beam.soll,dWf=w['dWf'])  # Trace3D
        section.add_element(cav)
    Phys['RFSection']=section.length
    return section  
def make_half_cell(w,upstream=True):  # 1/2 cell
    tki  = Beam.soll.tkin  
    gaps = w['gaps']
    ld   = w['ld']
    lcav = w['lcav']
    ld   = 0.5*(ld-2.*gaps*lcav)
    lqf  = 0.5*w['lqf']
    lqd  = 0.5*w['lqd']
    kq   = k0(gradient=w['dBdz'],tkin=tki)                  # quad strength @ entrance
    mqf=QF(k0=kq,length=lqf,label='QF',beam=Beam.soll)      # F quad 
    mqd=QD(k0=kq,length=lqd,label='QD',beam=Beam.soll)      # D quad 
    md=D(ld,beam=Beam.soll)                                 # drift zw. F quad und cavities
    cell=Lattice()
    if upstream : # 1/2 basis zelle upstream
        cell.add_element(mqf)   # QF
        cell.add_element(md)    # D
        rf_section = make_rf_section(w) # cavities
        cell.append(rf_section) # RF
        # -- energy update --
        mqd = mqd.update()
        md  = md.update()
        cell.add_element(md)    # D
        cell.add_element(mqd)   # QD
    else:  #  1/2 basis zelle downstream  = reverse of upstream
        cell.add_element(mqd)   # QD
        cell.add_element(md)    # D
        rf_section = make_rf_section(w) # cavities
        cell.append(rf_section) # RF
        # -- energy update --
        md  = md.update()
        mqf = mqf.update()
        cell.add_element(md)    # D
        cell.add_element(mqf)   # QF
    Phys['LQF']=2.*mqf.length
    Phys['LQD']=2.*mqd.length
    Phys['LD'] =md.length
    Phys['CELL']=cell.length
    deltaTK=Beam.soll.tkin - tki
    return cell,deltaTK
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
def loesung():  # total classic FODO lattice
    # längen
    lqd  =  0.4     # QD len
    lqf  =  0.4     # QF len
    ld   =  0.4     # drift len                # KNOB effective focus of FODO
    # cavity werte
    lcav   =  Phys['spalt_laenge']         
    u0     = Phys['spalt_spannung']
    phi0   = Phys['soll_phase']*Phys['radians']
    fRF0   = Phys['frequenz']
    gaps_per_half_cell= 3                      # KNOB  gaps/cell
    dWf=1.                                     # acceleration flag
    # beam werte
    Beam.soll = Beam(Phys['injection_energy'])
    Phys['sigx_i'] = 5.e-3                     # KNOB sigma x (i)
    Phys['sigy_i'] = 2.5e-3                    # KNOB sigma y (i)
    Phys['dP/P']   = 2.e-2                     # KNOB dp/p (i)
    tk0 = Beam.soll.tkin                       # KNOB injection energy
    # fokusierung
    # dBdz0  = Phys['quad_gradient']*7.85      # KNOB quad gradient
    dBdz0  = Phys['quad_gradient']*8.2         # KNOB quad gradient
    # struktur werte
    ring = False                                # KNOB ring or transfer ?
    nboff_super_cells = 16*10                  # KNOB  final energy
    # nboff_super_cells = 16*5                 # KNOB  final energy
    nboff_super_cells = 16*1                   # KNOB  final energy
    # nboff_super_cells = 3                    # KNOB  final energy
    w ={'lqd':lqd,
        'lqf':lqf,
        'ld':ld,
        'lcav':lcav,
        'U0':u0,
        'phi0':phi0,
        'fRF':fRF0,
        'dBdz':dBdz0,
        'dWf':dWf,
        'gaps':gaps_per_half_cell,}
    #--------------------1st cavity-----------
    gapi = RFG(U0=w['U0'],PhiSoll=w['phi0'],fRF=w['fRF'],beam=Beam.soll)  # Trace3D
    s_ttf_i =gapi.tr
    # objprnt(gapi, filter={'matrix'})
    #-----------------------------------------
    super_cell=Lattice()
    nboff_gaps=0                 # gap counter
    for icell in range(nboff_super_cells):
        zelle = Lattice()  # basis zelle
        (half_cell,deltaTK) = make_half_cell(w,upstream=True); nboff_gaps+=gaps_per_half_cell
        zelle.append(half_cell)
        (half_cell,deltaTK) = make_half_cell(w,upstream=False); nboff_gaps+=gaps_per_half_cell
        zelle.append(half_cell)
        super_cell.append(zelle)  # add zelle to super cell
    # super_cell.out()
    lattice_length=super_cell.length
    # print('lattice length [m]={}'.format(lattice_length))
    #-----------------------------------------
    # Berechne ganze Zelle und Anfangswerte 
    mcell,betax,betay = super_cell.cell(closed=ring)
    print('\nBETAx(i) {:.3g} [m], BETAy(i) {:.3g} [m]'.format(betax,betay))
    #---------------last cavity---------------
    gapf = RFG(U0=w['U0'],PhiSoll=w['phi0'],fRF=w['fRF'],beam=Beam.soll)  # Trace3D
    s_ttf_f =gapf.tr
    # objprnt(gapf, filter={'matrix'})    
    #-----------------------------------------
    # Zusammenfassung
    s_tk_i  =tk0
    s_tk_f  =Beam.soll.tkin
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
    s_emi   =(Phys['emitx_i'],Phys['emity_i'])
    s_aper  =(Phys['sigx_i'],Phys['sigy_i'])
    s_phis  =Phys['soll_phase']
    s_lamb  =Phys['wellenlänge']
    s_freq  =Phys['frequenz']
    s_nboff_cells =nboff_super_cells
    summary={
    'quadrupole size          [m]':s_lqd,
    'particle rest mass[MeV/c**2]':s_e0,
    'particle energy(i,f)   [Mev]':(s_tk_i,s_tk_f),
    'quad strength(i,f)  [1/m**2]':(s_kq_i,s_kq_f),
    'number of cavities          ':s_gaps,
    'qudrupole gradient     [T/m]':s_bgrad,
    'av. voltage_per_gap     [MV]':s_u0,
    'lattice_length           [m]':s_latlen,
    'tot.acceleration       [MeV]':s_utot,
    'av.acceleration      [MeV/m]':s_accel,
    'emittance-x,-y(i)    [m*rad]':s_emi,
    'sigma-x,-y(i)            [m]':s_aper,
    'sync. phase            [deg]':s_phis,
    'wave length              [m]':s_lamb,
    'frequency              [MHz]':s_freq,
    'time transition factor (i,f)':(s_ttf_i,s_ttf_f),
    'number of cells             ': s_nboff_cells,
    }
    dictprnt(summary,'Summary')
    #-----------------------------------------
    # Grafik: lösungen als f(s)
    functions = super_cell.functions(30)   
    display(functions)
if __name__ == '__main__':
    loesung()