import os, sys, math
"""
Evaluate various formulas from T.Wangler's book (TW)
"""
frq   = 800.e6    # [Hz] frequency
cl    = 2.998e8   # [m/s] speed of light
m0c2  = 938.27e6  # [eV] proton rest mass
tkin  = 50.       # [MeV] kin. energy
gap   = 0.048     # [m] rf-gap length
rbore = 0.020     # [m] rf-gap bore radius

pi    = math.pi
twopi = 2.*pi
omega = twopi*frq # [rad/s] phase frequency
W     = m0c2 + tkin*1.e6  # [eV]
gamma = W/m0c2
beta  = math.sqrt(1.-1./gamma**2)
lamb  = cl/frq    # [m]

def I0(x):
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
        # DEBUG_MODULE('(I0,x )',(res,x))
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
            # DEBUG_MODULE('(I0,x )',(res,x))
        except OverflowError as ex:
            print('Bessel-function I0 overflow: (arg = {:6.3f})! - STOP'.format(x))
            sys.exit(1)
    return res

def TWpp28():
    Rc =     2.405*cl/omega     # TW pp. 28
    print('TW-pp.28: frequency(frq)[Mhz]= {:5.3f}, pillbox radius(Rc)[cm]= {:5.3f}'.format(frq*1.e-6,Rc*1.e2))

def TW2_24():
    K2    = (twopi/lamb)**2*(1.-1./beta**2)   # TW 2.24
    K     = math.sqrt(abs(K2))
    print('TW-2.24 : K2= {}, K= {}'.format(K2,K))

def TW2_43():
    K2    = (twopi/lamb)**2*(1.-1./beta**2)   # TW 2.24
    K     = math.sqrt(abs(K2))
    x     = (pi*gap)/(beta*lamb)
    Ttfk  = 1./I0(K*rbore)*math.sin(x)/x   # TW 2.43
    print('TW-2.43 : T(k)~= {}'.format(Ttfk))
    
if __name__ == '__main__':
    TWpp28()
    TW2_24()
    TW2_43()