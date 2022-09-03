#!/Users/klotz/SIMULINAC_env/bin/python
# -*- coding: utf-8 -*-
"""
Copyright 2015 Wolf-Dieter Klotz <wdklotz@gmail.com>
This file is part of the SIMULINAC code

    SIMULINAC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    SIMULINAC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SIMULINAC.  If not, see <http://www.gnu.org/licenses/>.
"""
import os, sys, math
"""
Evaluate various formulas from books/articles: 
T.Wangler (TW), 
H.Wiedemann (HW), 
E.Jensen (EJ) CERN

"""
frq     = 816.e6    # [Hz] frequency
cl      = 2.998e8   # [m/s] speed of light
m0c2    = 938.27e6  # [eV] proton rest mass
tkin    = 50.       # [MeV] kin. energy
gap     = 0.048     # [m] rf-gap length
rbore   = 0.020     # [m] rf-gap bore radius
inv_sig = 1.7e-8    # [ohms*m]
E0      = 5.e6      # [MV/m] peak cavity field

pi    = math.pi
twopi = 2.*pi
omega = twopi*frq # [rad/s] phase frequency
W     = m0c2 + tkin*1.e6  # [eV]
gamma = W/m0c2
beta  = math.sqrt(1.-1./gamma**2)
lamb  = cl/frq    # [m]

def I0(x):      #TODO   kw
    """
    Modified Bessel function I of integer order 0
    ref.: Hanbook of Mathematical Functions, M.Abramowitz & I.A.Stegun
    """
    t = x/3.75
    if 0. <= x and x < 3.75:
        t2 = t*t
        res = 1.
        res+= 3.5156229*t2
        res+= 3.0899424*t2*t2
        res+= 1.2067492*t2*t2*t2
        res+= 0.2659732*t2*t2*t2*t2
        res+= 0.0360768*t2*t2*t2*t2*t2
        res+= 0.0045813*t2*t2*t2*t2*t2*t2
        print('(I0,x )',(res,x))
    elif 3.75 <= x:
        tm1 = 1./t
        res = 0.39894228
        res+= 0.01328529*tm1
        res+= 0.00225319*tm1*tm1
        res-= 0.00157565*tm1*tm1*tm1
        res+= 0.00916281*tm1*tm1*tm1*tm1
        res-= 0.02057706*tm1*tm1*tm1*tm1*tm1
        res+= 0.02635537*tm1*tm1*tm1*tm1*tm1*tm1
        res-= 0.01647633*tm1*tm1*tm1*tm1*tm1*tm1*tm1
        res+= 0.00392377*tm1*tm1*tm1*tm1*tm1*tm1*tm1*tm1
        try:
            res = res*exp(x)/sqrt(x)
            print('(I0,x )',(res,x))
        except OverflowError as ex:
            print('Bessel-function I0 overflow: (arg = {:6.3f})! - STOP'.format(x))
            sys.exit(1)
    return res
def J1(x):
    """ -3 <= x <= +3 """
    xov3 = x/3.
    j1  = 0.5 
    j1 -= 0.56249985*xov3**2
    j1 += 0.21093573*xov3**4
    j1 -= 0.03954289*xov3**6
    j1 += 0.00443319*xov3**8
    j1 -= 0.00031761*xov3**10
    j1 += 0.00001109*xov3**12
    j1  = x*j1
    return j1
def TWpp28():   # TW pp. 28
    Rc =     2.405*cl/omega
    print('TW-pp.28 : frequency(frq)[Mhz]= {:5.3f}, pillbox radius(Rc)[cm]= {:5.3f}'.format(frq*1.e-6,Rc*1.e2))
    return Rc
def TW2_24():   # TW 2.24
    K2    = (twopi/lamb)**2*(1.-1./beta**2)
    K     = math.sqrt(abs(K2))
    print('TW-2.24  : K2= {}, K= {}'.format(K2,K))
def TW2_43():   # TW 2.43
    K2    = (twopi/lamb)**2*(1.-1./beta**2)
    K     = math.sqrt(abs(K2))
    x     = (pi*gap)/(beta*lamb)
    Ttfk  = 1./I0(K*rbore)*math.sin(x)/x
    print('TW-2.43  : T(k)~= {} - time transition factor'.format(Ttfk))
    return Ttfk
def HWII6_61(): # HW vol II 6.61
    rs = 1.28 * math.sqrt(frq)
    print('HWII-6.61: rs[MOhm/m]~= {} - specific shunt impedance (pillbox)'.format(rs))
def EJ_26():   # EJ 26
    eta   = 376.73    # Ohms
    chi01 = 2.40483
    fact  = 4.*eta/(chi01**3*pi*J1(chi01)**2)
    RA    = 8.e-3     # [Ohm] Cu surface resistance @ 1 GHz scales with sqrt(omega)
    sig   = 5.8e7     # [S/m] Cu conductivity
    E0LT  = E0*gap*TW2_43()
    Vacc  = E0LT
    # Alceli pillbox
    l     = 44.e-3
    a     = TWpp28()
    lOa   = l/a 
    ROQ   = fact*math.sin(chi01/2*lOa)**2/lOa   # R/Q
    Q0    = math.sqrt(2.*a*eta*sig*chi01)/(2.*(1.+1./lOa))  # Q  EJ 37
    P     = Vacc**2/2/(ROQ*Q0)  # [W] power EJ 32
    print('EJ-26    : R/Q[Ohm]= {:4.4g} Q= {:1.4e} (Alceli pillbox l/a ={})'.format(ROQ,Q0,lOa))
    print('         : Vacc[KV]= {:4.4g} P[KW]= {:4.4g}'.format(Vacc*1.e-3,P*1.e-3))
    # DESY 500MHz pillbox
    dsyRC = 462./2.*1e-3   # [m]
    dsyL  = 276.*1e-3      # [m]
    l    = dsyL
    a    = dsyRC
    lOa  = l/a
    ROQ  = fact*math.sin(chi01/2*lOa)**2/lOa   # R/Q
    Q0   = math.sqrt(2.*a*eta*sig*chi01)/(2.*(1.+1./lOa))  # Q  EJ 37
    P    = Vacc**2/2/(ROQ*Q0)  # [W] power
    print('         : R/Q[Ohm]= {:4.4g} Q= {:1.4e} (desy pillbox l/a={})'.format(ROQ,Q0,lOa))
    print('         : Vacc[KV]= {:4.4g} P[KW]= {:4.4g}'.format(Vacc*1.e-3,P*1.e-3))
    
if __name__ == '__main__':
    TW2_24()
    TW2_43()
    HWII6_61()
    EJ_26()