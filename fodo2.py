#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import CONF,k0,dictprnt,objprnt,Beam,Proton,Electron
from elements import D,QF,QD,RFG,RFC
from lattice import Lattice
from pylab import plot,show,legend,figure,subplot,axis
from math import sqrt,radians
import fileLoader

def display(functions):                 ## plotting
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
    bxn = [-x for x in bx]    # beta-x (negatif)
    byn = [-x for x in by]    # beta-y (negatif)
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
    viseox = [x*vscale for x in viseo]
    plot(z,viseox,label='',color='black')
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
def make_rf_section(w):                 ## RF sektion
    gaps = w['gaps']    # gaps/half-cell
    section = Lattice()
    for i in range(gaps):
        cav=RFC(length=w['lcav'],U0=w['U0'],PhiSoll=w['phi0'],fRF=w['fRF'],beam=Beam.soll,dWf=w['dWf'])  # Trace3D
        # objprnt(Beam.soll)
        section.add_element(cav)
    CONF['RFSection']=section.length
    return section  
def make_half_cell(w,upstream=True):    ## 1/2 cell
    tki  = Beam.soll.tkin  
    gaps = w['gaps']
    ld   = w['ld']
    lcav = w['lcav']
    ld   = 0.5*(ld-gaps*lcav)
    if ld < 0.:
        raise RuntimeWarning('negative drift space!')
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
    CONF['LQF']=2.*mqf.length
    CONF['LQD']=2.*mqd.length
    CONF['LD'] =md.length
    CONF['CELL']=cell.length
    deltaTK=Beam.soll.tkin - tki
    return cell,deltaTK
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
def loesung():                          ## total classic FODO lattice (1st result, used as reference!)
    ring = True                                # KNOB ring or transfer ?
    #-----------------------------------------
    print('Injected beam:\n'+Beam.soll.out(False))
    super_cell = fileLoader.read_yaml_and_parse()
    #-----------------------------------------
    # Berechne ganze Zelle und Anfangswerte 
    mcell,betax,betay = super_cell.cell(closed=ring)
    print('\nBETAx(i) {:.3g} [m], BETAy(i) {:.3g} [m]'.format(betax,betay))
    #-----------------------------------------
    # Grafik: lösungen als f(s)
    functions = super_cell.functions(30)   
    display(functions)
if __name__ == '__main__':
    loesung()