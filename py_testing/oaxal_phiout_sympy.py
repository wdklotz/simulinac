from sympy import *

print("======================Uebung======================")
a,b,c,d,e, f = symbols('a b c d e f')
x, y = symbols('x y')

p=a+b*x+c*y+d*x*y+e*x**2+f*y**2

print('polynom= ',p,'\n')
ex1 = p.subs(x,0)
print('nur y-terme ',ex1)  
print('absoluter term ',ex1.subs(y,0))
print('linearer y-term ',ex1.coeff(y,1),'\n') 

print("==============================================================================")
print("Use SYMPY to make the error-prone factorization work for the Open-XAL formulas")
print("==============================================================================")
c         = symbols('c')
mc02      = symbols('mc02')
m0c3      = symbols('m0c3')
gammas    = symbols('gammas')
betas     = symbols('betas')
qV0       = symbols('qV0')        # qE0L
omega     = symbols('omega')
db2bs     = symbols('db2bs')
DPHI      = symbols('DPHI')
DPHIS     = symbols('DPHIS')
dphi      = symbols('dphi')
phis      = symbols('phis')

cphis, sphis = symbols('cphis sphis')
cphi, sphi   = symbols('cphi sphi')
Tks, Tpks, Tppks   = symbols('Tks Tpks Tppks')   # T(ks), T'(ks) fuer SOLL
Sks, Spks, Sppks   = symbols('Sks Spks Sppks')   # S(ks), S'(ks) fuer SOLL
Tk, Tpk, Tppk      = symbols('Tk Tpk Tppk')   # T(ks), T'(ks) fuer SOLL
Sk, Spk, Sppk      = symbols('Sk Spk Sppk')   # S(ks), S'(ks) fuer SOLL
FAC   = symbols('FAC')
FAC1  = symbols('FAC1')
FAC2  = symbols('FAC2')
_tpk  = symbols('_tpk')
_spk  = symbols('_spk')
_cphi = symbols('_cphi')
_sphi = symbols('_sphi')
_fac  = symbols('_fac')
_fac1 =symbols('_fac1')
gb3   = symbols('gb3')
_fac2 = symbols('_fac2')
_gb3  = symbols('_gb3')
gsbs3 = symbols('gsbs3')


init_printing()

DPHIS = FAC1*(Tpks*sphis+Spks*cphis)
print('DPHIS= \n',DPHIS,'\n')
_fac1 = qV0*omega/m0c3*gsbs3
DPHIS=DPHIS.subs(FAC1,_fac1)
print('DPHIS= \n',DPHIS,'\n')

print('DPHIS= \n',expand(DPHIS),'\n')

DPHI = FAC*(Tpk*sphi+Spk*cphi)
print('DPHI= \n',DPHI,'\n')

_tpk=Tpks-Tppks*FAC2*db2bs
_spk=Spks-Sppks*FAC2*db2bs
DPHI = DPHI.subs(Tpk,_tpk)
DPHI = DPHI.subs(Spk,_spk)
print('DPHI= \n',DPHI,'\n')
_cphi=cphis-sphis*dphi
_sphi=sphis+cphis*dphi
DPHI=DPHI.subs(cphi,_cphi)
DPHI=DPHI.subs(sphi,_sphi)
print('DPHI= \n',DPHI,'\n')
_fac=qV0*omega/m0c3*gb3
DPHI=DPHI.subs(FAC,_fac)
print('DPHI= \n',DPHI,'\n')
_fac2 = omega/c/betas*db2bs
DPHI=DPHI.subs(FAC2,_fac2)
print('DPHI= \n',DPHI,'\n')
_gb3 = gsbs3-3/gammas/betas**3*db2bs
DPHI=DPHI.subs(gb3,_gb3)
print('DPHI= \n',DPHI,'\n')

print('DPHI= \n',expand(DPHI),'\n')

DDPHI = symbols('DDPHI')
DDPHI = DPHI - DPHIS
print('DDPHI= \n',DDPHI,'\n')

print('DDPHI(expanded)= \n',expand(DDPHI),'\n')

