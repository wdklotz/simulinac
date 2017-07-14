from math import sqrt

import elements as elm
from lattice import Lattice
from setutil import Proton, DEBUG, dictprnt, dBdxprot

def m_fodo(params):
    kq       = params['kq [1/m^2]']
    lq       = params['quad [m]']
    ld       = params['drift [m]']
    particle = params['particle']
    qf = elm.QF(k0=kq, length=lq/2., label='QF/2', particle=particle)
    qd = elm.QD(k0=kq, length=lq,    label='QD',   particle=particle)
    dr = elm.D(        length=ld,    label='D',    particle=particle)
    
    cell = Lattice()
    cell.add_element(qf)
    cell.add_element(dr)
    cell.add_element(qd)
    cell.add_element(dr)
    cell.add_element(qf)
    return cell
    
def test0(params):
    print('------------------------------Test0--')
    ld       = params['drift [m]']
    lq       = params['quad [m]']     # full quad
    tkin     = params['tk [Mev]']
    grad     = params['dB/dx [T/m]']

    L        = ld+2.*lq               # 1/2 cell
    particle = Proton(tkin)
    kq       = grad/particle.brho
    f        = 1./(kq*lq) # focus [m]
    kappa    = f/L   # >1

    params['kq [1/m^2]']= kq
    params['B*rho [Tm]']= particle.brho
    params['particle']  = particle
    params['focus [m]'] = f
    params['kappa']     = kappa
    
    cell = m_fodo(params)
    params['cell [m]']  = cell.length
    cell.cell()    
    dictprnt(params,'params from input')
    # DEBUG('cell:',cell.string())
    
    # params with optimal kappa (Wiedemann I,pp189,(6.16)
    kappa = sqrt(2.)
    f = kappa*L
    kq = 1./(f*lq)
    grad = dBdxprot(kq,tkin)
    params['dB/dx [T/m]'] = grad
    params['kq [1/m^2]']  = kq
    params['B*rho [Tm]']  = particle.brho
    params['focus [m]']   = f
    params['kappa']       = kappa
    
    cell = m_fodo(params)
    cell.cell()    
    dictprnt(params,'params with optimal kappa (Wiedemann I,pp189,(6.16)')

#-------------------------------------------main---
if __name__ == '__main__':
    params = {
        'drift [m]'  : 0.6,
        'quad [m]'   : 0.1,
        'tk [Mev]'   : 5.0,
        'dB/dx [T/m]': 9.72
        }
    test0(params)
    