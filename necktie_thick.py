import elements as elm
from lattice import Lattice
from setutil import Proton, DEBUG, dictprnt

def test0(params):
    print('------------------------------Test0--')
    ld       = params['drift [m]']
    lq       = params['quad [m]']     #full quad
    tkin     = params['tk [Mev]']
    grad     = params['dB/dx [T/m]']
    particle = Proton(tkin)
    kq       = grad/particle.brho

    qf = elm.QF(k0=kq, length=lq/2., label='QF/2', particle=particle)
    qd = elm.QD(k0=kq, length=lq,    label='QD',   particle=particle)
    dr = elm.D(        length=ld,    label='D',    particle=particle)
    
    cell = Lattice()
    cell.add_element(qf)
    cell.add_element(dr)
    cell.add_element(qd)
    cell.add_element(dr)
    cell.add_element(qf)

    params['kq [1/m^2]']= kq
    params['B*rho [Tm]']= particle.brho
    params['particle']  = particle.name
    params['cell [m]']  = cell.length
    
    dictprnt(params,'params')
    DEBUG('cell:',cell.string())

    cell.cell()
#-------------------------------------------main---
if __name__ == '__main__':
    params = {
        'drift [m]'  : 0.6,
        'quad [m]'   : 0.1,
        'tk [Mev]'   : 5.,
        'dB/dx [T/m]': 9.72
        }
    test0(params)
    