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

qV0, omega, c, betas, db2bs, dphi = symbols('qV0 omega c betas db2bs dphi')
Tks, Sks, Tpks, Spks, cphis, sphis  = symbols('Tks Sks Tpks Spks cphis sphis')

_TK=Tks - Tpks * omega/c/betas*db2bs
_SK=Sks - Spks * omega/c/betas*db2bs
_CPHI=cphis-sphis*dphi
_SPHI=sphis+cphis*dphi
_dw  = qV0*(_TK*_CPHI-_SK*_SPHI)
_dws = qV0*(Tks*cphis-Sks*sphis)
print('dw= \n',_dw,'\n')
print('dw(expanded)= \n',expand(_dw),'\n')
print('dws= \n',_dws,'\n')

Wins, dws, Win, dw, Wout, Wouts, DWout = symbols('Wins dws Win dw Wout Wouts DWout')
DWout2Wout = DWout/Wout
Wout       = Win+dw
Wouts      = Wins+dws
print('Wouts= \n',Wouts,'\n')
print('Wout= \n',Wout,'\n')
# print('DWout/Wout= \n',DWout2Wout,'\n')

# _Wouts = Wins+dws
# _Wout  = Win+dw
DWout = Wout-Wouts
DWout = DWout.subs(dw,_dw)
DWout = DWout.subs(dws,_dws)
# DWout = DWout.subs(dw,_dw)
# DWout = DWout.subs(dws,_dws)
print('DWout= \n',DWout,'\n')

DWout=expand(DWout)
print('DWout(expanded)= \n',DWout,'\n')

# delta_beta_terms = DWout.subs(dphi,0)
# print('delta_beta_terms= \n',delta_beta_terms,'\n')

# absolute_terms = delta_beta_terms.subs(db2bs,0)
# print('absolute_terms= \n',absolute_terms,'\n')

# linear_delta_beta_terms = delta_beta_terms.coeff(db2bs,1)
# print('linear_delta_beta_terms= \n',linear_delta_beta_terms,'\n')

# delta_phi_terms = DWout.subs(db2bs,0)
# print('delta_phi_terms= \n',delta_phi_terms,'\n')

# linear_delta_phi_terms = delta_phi_terms.coeff(dphi,1)
# print('linear_delta_phi_terms= \n',linear_delta_phi_terms,'\n')

# print('linear_delta_beta_terms(simplified)= \n',simplify(linear_delta_beta_terms),'\n')
# print('linear_delta_phi_terms(simplified)= \n',simplify(linear_delta_phi_terms),'\n')

gamma, m0c2 = symbols('gamma m0c2')
Dp2pout = gamma/(gamma+1)*DWout/Wout
_gamma = 1.+ Wout/m0c2
Dp2pout = Dp2pout.subs(gamma,_gamma)
print('dP/p out= \n',simplify(Dp2pout))
