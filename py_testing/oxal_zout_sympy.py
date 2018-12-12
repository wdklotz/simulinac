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
# mc02       = symbols('mc02')
m0c3       = symbols('m0c3')
gammas       = symbols('gammas')
betas      = symbols('betas')
# beta       = symbols('beta')
db2bs     = symbols('db2bs')
dphi      = symbols('dphi')
qV0       = symbols('qV0')        # qE0L
omega       = symbols('omega')
# z        = symbols('z')
phin      = symbols('phin')
phis      = symbols('phis')

cphis, CPHI  = symbols('cphis CPHI')      # cos(phis) fuer SOLL
sphis, SPHI  = symbols('sphis SPHI')      # sin(phis) fuer Soll

fac, GBHM3, gbshm3 = symbols('fac GBHM3 gbshm3')
Tks, Tpks, Tppks   = symbols('Tks Tpks Tppks')   # T(ks), T'(ks) fuer SOLL
Tk, TPK, Tppk      = symbols('Tk  TPK  Tppk')   # T(ks), T'(ks) fuer SOLL
Sks, Spks, Sppks   = symbols('Sks Spks Sppks')   # S(ks), S'(ks) fuer SOLL
Sk, SPK, Sppk      = symbols('Sk  SPK  Sppk')   # S(ks), S'(ks) fuer SOLL

init_printing()

BET  = symbols('BET')
DPHI = symbols('DPHI')
DPHIS = symbols('DPHIS')
zout = -c*BET/omega*(phin + DPHI - (phis+DPHIS))
# print('zout= ',zout)

FAC   = symbols('FAC')
# DPHI = FAC*(TPK*SPHI-SPK*CPHI)
zout=zout.subs(DPHI,FAC*(TPK*SPHI+SPK*CPHI))
zout=zout.subs(DPHIS,FAC*(Tpks*sphis+Spks*cphis))
print('zout= \n',zout,'\n')

# CPHI  = cphis-sphis*dphi      # cos(phi) = cos(phis)-sin(phis)*deltaPhi
# SPHI  = sphis+cphis*dphi      # selbe Naeherung wie fuer cos
zout = zout.subs(CPHI,cphis-sphis*dphi)
zout = zout.subs(SPHI,sphis+cphis*dphi)
# print('zout= ',zout)

# TPK   = Tpks - Tppks*db2bs*omg/(c*bets)
# SPK   = Spks - Sppks*db2bs*omg/(c*bets)
zout = zout.subs(TPK,Tpks - Tppks*db2bs*omega/(c*betas))
zout = zout.subs(SPK,Spks - Sppks*db2bs*omega/(c*betas))
# print('zout= ',zout)

# FAC = qV0*omg/mc3*GBHM3
zout = zout.subs(FAC,qV0*omega/mc03*GBHM3)
# print('zout= ',zout)

# BET = bets*(1+db2bs)
zout = zout.subs(BET,betas*(1+db2bs))
# print('zout= ',zout)

# GBHM3 = gbshm3-3/gas/bets**3*db2bs
zout = zout.subs(GBHM3,gbshm3-3/gammas/betas**3*db2bs)
# print('zout= ',zout,'\n')

zout= expand(zout)
print('expand(zout)= \n',zout,'\n')

# delta_beta_terms = zout.subs(dphi,0)
# print('DeltaBeta terms ',delta_beta_terms,'\n')

# absolute_terms_1 = delta_beta_terms.subs(db2bs,0)
# print('Absolute terms= \n',simplify(absolute_terms_1),'\n')

# linear_deltabeta_terms = delta_beta_terms.coeff(db2bs,1)
# print('Linear DeltaBeta terms= \n ',linear_deltabeta_terms,'\n')

# delta_phi_terms = zout.subs(db2bs,0)
# print('DeltaPhi terms ',delta_phi_terms,'\n')

# absolute_terms_0 = delta_phi_terms.subs(dphi,0)
# print('Absolute terms= \n ',absolute_terms_0,'\n')

# linear_deltaphi_terms = delta_phi_terms.coeff(dphi,1)
# print('Linear DeltaPhi terms= \n ',linear_deltaphi_terms,'\n')

# print('Absolute terms(simplified)= \n ',simplify(absolute_terms_0),'\n')
# print('Delta_beta_terms(simplified)= \n ',simplify(linear_deltabeta_terms),'\n')
# print('Delta_phi_terms(simplified)= \n ',simplify(linear_deltaphi_terms),'\n')

