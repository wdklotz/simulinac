#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import Phys
from elements import k0,I,D,QF,QD,SD,WD,CAV,RFG,Beam
from lattice import Lattice
from pylab import plot, show, legend, figure, subplot, axis
from math import sqrt

def display(functions,title):  ## plotting
    #----------*----------*   # unpack
    beta_fun = functions[0]
    cos_like = functions[1]
    sin_like = functions[2]
    emitx=Phys['emitx(i)']  # emittance @ entrance
    emity=Phys['emity(i)']  # emittance @ entrance
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
def make_cavity():   
    # cavity
    lcav   = Phys['spalt_laenge']
    u0     = Phys['spalt_spannung']*0.05
    phi0   = Phys['soll_phase']*Phys['radians']
    fRF0   = Phys['frequenz']
    dWf    = 1.

    tk = Beam.soll.tkin                    # kinetic energy @ entrance
    l=0.5*lcav
    cavity = Lattice()
    dri = D(length=0.5*l,beam=Beam.soll,label='>')   # drift before RFgap
    gap = RFG(U0=u0,PhiSoll=phi0,fRF=fRF0,label='rfg',beam=Beam.soll,dWf=dWf)  # Trace3D
    drf = D(length=0.5*l,beam=Beam.soll,label='<')   # drift after RFgap
    cavity.add_element(dri)
    cavity.add_element(gap)
    cavity.add_element(drf)
    Phys['LCAV']=cavity.length
    # print('cavity_length: ',cavity.length)
    # cavity.out()
    # objprnt(Beam.soll,'gap exit',{'matrix'})
    return cavity
def make_rf_section(gaps=1):   ## RF sektion
    ''' gaps: nboff gaps per rf section'''
    section = Lattice()
    for i in range(gaps):
        cav = make_cavity()
        section.append(cav)
    Phys['RFSection']=section.length
    print('rf_section_length: ',section.length)
    return section  
def objprnt(what,text='========',filter={}):  ## helper to print objects as dictionary
        print('\n========= '+text+' =================')
        for k,v in what.__dict__.items():
            if k in filter:
                continue
            print(k.rjust(30),':',v)
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
    # k werte
    sk=0.6
    # sk=x
    kf=7.5*sk
    kd=7.0625*sk
    # längen
    lq=0.2
    ld=0.6
    # ld=x
    # rf
    gaps = 6
    # tk0    = Phys['kinetic_energy']*1.       # KNOB: injection energy

    rf_section = make_rf_section(gaps=gaps)
    rf_len = rf_section.length
    lhalf=0.5*ld
    l1=lhalf-rf_len
    print('rf_len {}, lhalf {}, l1 {}'.format(rf_len,lhalf,l1))
    
    mD1=D(length=l1,label='DR')
    
    lat=Lattice()
    lat.append(rf_section)   # RF
    lat.add_element(mD1)     # D
    lat.add_element(mQF1)    # F1
    lat.add_element(mQD1)    # D1    
    lat.add_element(mQF1)    # F1
    lat.add_element(mD1)     # D
    lat.append(rf_section)   # RF
    lat.append(rf_section)   # RF
    lat.add_element(mD1)     # D
    lat.add_element(mQD2)    # D2
    lat.add_element(mQF2)    # F2
    lat.add_element(mQD2)    # D2
    lat.add_element(mD1)     # D
    lat.append(rf_section)   # RF
    lat.append(rf_section)   # RF
    lat.add_element(mD1)     # D
    lat.add_element(mQF1)    # F1
    lat.add_element(mQD1)    # D1   
    lat.add_element(mQF1)    # F1
    lat.add_element(mD1)     # D
    lat.append(rf_section)   # RF
    
    lat.append(lat) # 2 cells
    lat.append(lat) # 4 cells
    lat.append(lat) # 8 cells
    
    mcell,betax,betay=lat.cell(closed=True)
    # mcell.out()
    print('energy(f)= {} [MeV]'.format(Beam.soll.tkin))
    
    Phys['sigx(i)']=sqrt(Phys['emitx(i)']*betax)
    Phys['sigy(i)']=sqrt(Phys['emity(i)']*betay)
    functions = lat.functions(40)   
    display(functions,'x {}'.format(x))
    return
#----------*----------*----------*----------*----------*----------*
if __name__ == '__main__':
    test0()
    # for x in [0.0+n*0.1 for n in range(30)]:
        # test0(x)
