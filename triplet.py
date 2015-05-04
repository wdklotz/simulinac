#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import Phys,dictprnt,objprnt
from elements import k0,I,D,QF,QD,SD,WD,CAV,RFG,Beam
from lattice import Lattice
from pylab import plot, show, legend, figure, subplot, axis
from math import sqrt

def display(functions,title):  ## plotting
    #----------*----------*   # unpack
    beta_fun = functions[0]
    cos_like = functions[1]
    sin_like = functions[2]
    emitx=Phys['emitx_i']  # emittance @ entrance
    emity=Phys['emity_i']  # emittance @ entrance
    #----------*----------*   # bahnkoordinate z
    z   = [ x[0] for x in beta_fun]    
    #----------*----------*
    bx  = [ sqrt(x[1]*emitx) for x in beta_fun]    # envelope (beta-x)
    by  = [ sqrt(x[2]*emity) for x in beta_fun]    # envelope (beta-y)
    bxn = [-x for x in bx]    # beta-x (negatif)
    byn = [-x for x in by]    # beta-y (negatif)
    zero= [0. for x in beta_fun]  # zero line
    #----------*----------*   # trajectories
    cx = [x[0] for x in cos_like]   # cos-like-x
    cy = [x[2] for x in cos_like]   # cos-like-y
    sx = [x[0] for x in sin_like]   # sin-like-x
    sy = [x[2] for x in sin_like]   # sin-like-x
    #----------*----------*   # figure frame
    figure(str(title))
    #----------*----------*   # transverse X
    splot=subplot(211)
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
    splot=subplot(212)
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
    #----------*----------*
    show(block=True)
def make_cavity(w):            ## one cavity
    # cavity
    lcav   = Phys['spalt_laenge']
    u0     = Phys['spalt_spannung']*0.05
    phi0   = Phys['soll_phase']*Phys['radians']
    fRF0   = Phys['frequenz']
    dWf    = w['dWf']
    cavity = Lattice()

    soll = Beam.soll                   # beam @ entrance
    tki = soll.tkin
    lcav05 = 0.5*lcav
    dri = D(length=lcav05, beam=soll, label='>')                               # drift before RFgap
    gap = RFG(U0=u0, PhiSoll=phi0, fRF=fRF0, label='rfg', beam=soll, dWf=dWf)  # Trace3D
    drf = D(length=lcav05, beam=soll, label='<')                               # drift after RFgap
    cavity.add_element(dri)
    cavity.add_element(gap)
    cavity.add_element(drf)
    # Phys['LCAV']=cavity.length
    # print('cavity_length: ',cavity.length)
    objprnt(soll,'beam @ gap exit',{'matrix'})
    return cavity
def make_rf_section(w):        ## many cavities
    gaps = w['gaps']
    section = Lattice()
    for i in range(gaps):
        cav = make_cavity(w)
        section.append(cav)
    # Phys['RFSection']=section.length
    print('rf_section_length: ',section.length)
    return section  
def make_cell(w):              ## cell
    gradient  = Phys['quad_gradient']
    cell      = Lattice()
    soll      = w['soll']
    ld        = w['ld']
    lq_short  = w['lq1']
    lq_long   = w['lq2']
    kf        = w['kf']
    kd        = w['kd']
    
    mQFs = QF(k0=kf,length=lq_short, label='QFs', beam=soll)
    mQFl = QF(k0=kf,length=lq_long,  label='QFl', beam=soll)
    mQDs = QD(k0=kd,length=lq_short, label='QDs', beam=soll)
    mQDl = QD(k0=kd,length=lq_long , label='QDl', beam=soll)
    mD   = D(length=ld, label='D', beam=soll)

    rf = make_rf_section(w)
    cell.append(rf)                       # RF
    
    # ---- update beam ENERGY --------
    mD   = mD.update()
    mQFs = mQFs.update()
    mQDl = mQDl.update()
    cell.add_element(mD)                  # D
    cell.add_element(mQFs)                # Fs
    cell.add_element(mQDl)                # Dl    
    cell.add_element(mQFs)                # Fs
    cell.add_element(mD)                  # D

    rf = make_rf_section(w)
    cell.append(rf)                       # RF
    rf = make_rf_section(w)
    cell.append(rf)                       # RF
    
    # ---- update beam ENERGY --------
    mD   = mD.update()
    mQFl = mQFl.update()
    mQDs = mQDs.update()
    cell.add_element(mD)                  # D
    cell.add_element(mQDs)                # Ds
    cell.add_element(mQFl)                # Fl
    cell.add_element(mQDs)                # Ds
    cell.add_element(mD)                  # D
    
    rf = make_rf_section(w)
    cell.append(rf)                       # RF
    return cell
def test0():
    # k werte
    sk=1.0
    # sk=x
    kf=7.5*sk
    kd=7.0625*sk
    # längen
    lq=0.2

    mQF1=QF(k0=kf,length=lq   ,label='F1')
    mQD1=QD(k0=kd,length=2.*lq,label='D1')
    mQF2=QF(k0=kd,length=2.*lq,label='F2')
    mQD2=QD(k0=kf,length=lq   ,label='D2')

    mFDF=mQF1*mQD1*mQF1 # 1st triplet
    mDFD=mQD2*mQF2*mQD2 # 2nd triplet
    mFDF.out()
    mDFD.out()
    fx1=-1./mFDF.matrix[1,0]
    fy1=-1./mFDF.matrix[3,2]
    fx2=-1./mDFD.matrix[1,0]
    fy2=-1./mDFD.matrix[3,2]
    print('focals triplet 1: fx= {:.3f}  fy= {:.3f}'.format(fx1,fy1))
    print('focals triplet 2: fx= {:.3f}  fy= {:.3f}'.format(fx2,fy2))
def test1(x):
    wert={}
    wert['dWf'] = 0.
    
    # beam
    tk0       = Phys['injection_energy']*1.       # KNOB: injection energy
    Beam.soll = Beam(tk0)
    wert['soll'] = Beam.soll
    
    # k werte
    sk=1.0
    # sk=x
    kf=7.5*sk
    kd=7.0625*sk
    wert['kf']=kf
    wert['kd']=kd

    # längen
    lq1 = 0.2
    lq2 = 2.*lq1
    ld=0.6
    # ld=x
    wert['lq1']=lq1
    wert['lq2']=lq2
    wert['ld']=ld

    # rf
    gaps = 3
    rf_len = gaps*Phys['spalt_laenge']
    lhalf=0.5*ld
    l1=lhalf-rf_len
    wert['ld']=l1
    wert['gaps'] = gaps

    # super cell
    ncells = 3
    super_cell = Lattice()

    dictprnt(Phys,'Phys')
    dictprnt(wert,'w')    
    
    for ncellsc in range(ncells):
        cell = make_cell(wert)
        super_cell.append(cell)
    
    mcell,betax,betay=super_cell.cell(closed=True)
    print('energy(f)= {} [MeV]'.format(Beam.soll.tkin))
    
    Phys['sigx_i']=sqrt(Phys['emitx_i']*betax)
    Phys['sigy_i']=sqrt(Phys['emity_i']*betay)
    functions = super_cell.functions(40)   
    display(functions,'x {}'.format(x))
    return
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    # test0()
    test1(0)
    # for x in [0.0+n*0.1 for n in range(30)]:
        # test0(x)
