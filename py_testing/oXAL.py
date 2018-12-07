from sympy import *

"""
Use SYMPY to make the error-prone factorization work for the Open-XAL formulas
"""
m0c2, gammas, betas, Ws, omega , c, Dp2p  \
          = symbols('m0c2 gammas betas Ws omega c Dp2p')
Tks, Tpks = symbols('Tks Tpks')   # T(ks), T'(ks) fuer SOLL
Sks, Spks = symbols('Sks Spks')   # S(ks), S'(ks) fuer SOLL
cphis     = symbols('cphis')      # cos(phis) fuer SOLL
sphis     = symbols('sphis')      # sin(phis) fuer Soll
qV0       = symbols('qV0')        # qE0L
omcb      = symbols('omcb')       # factor omega/c/betas
Dbeta     = symbols('Dbeta')      # deltaBeta = beta - beatas
Dphi      = symbols('Dphi')       # deltaPhi  = phi - phis
fac          = symbols('fac')        # factor
Tppks, Sppks = symbols('Tppks Sppks')   # 2-te Ableitungen T(k),S(k) nach k

init_printing()

cphi  = cphis-sphis*Dphi      # cos(phi) = cos(phis)-sin(phis)*deltaPhi
sphi  = sphis+cphis*Dphi      # selbe Naeherung wie fuer cos

Dk    = -omcb*Dbeta           # deltak = -omega/c/betas**2 * deltaBeta
Tk    = Tks-Tpks*omcb*Dbeta   # T(k) = T(ks) - T'(ks) * deltak
Sk    = Sks-Spks*omcb*Dbeta   # selbe Naeherung wie fuer T(k)
DWs   = qV0*(Tks*cphis-Sks*sphis)    # basis Formel (4.6.1) fuer SOLL
DW    = qV0*(Tk*cphi-Sk*sphi)        # basis Formel (4.6.1) fuer Particle
DWout    = DW-DWs                    # D(DW)  = DW(PARTICLE) - DW(SOLL)

print("-------------------delta W")
print(expand(DWout))

print("------------------delta beta")
Db1 = betas**(-1/2)*m0c2**(-1)*gammas**(-3)   # Konversion: deltaBeta <== deltaW
dw1 = gammas**(-1)*(gammas+1)*Ws*Dp2p         # Konversion: deltaW <== deltaP/P
Dbeta  = Db1*dw1                            # Konversion: deltaBeta(PARTICLE) <== deltaP/P
print(Dbeta)

print("------------delta phi")
Tpk = Tpks - Tppks * omcb * Dbeta    # (4.6.4)
Spk = Spks - Sppks * omcb * Dbeta    # (4.6.4)
DPs = fac*(Tpks*sphis+Spks*cphis)    # delta Phase SOLL   (4.6.2)
Dp  = fac*(Tpk*sphi+Spk*cphi)      # delta Phase PARTICLE (4.6.2)
DPout = Dp-DPs
print(expand(DPout))

print("------------------ombc=omega/c/bstas**2")
ombc = omega/c/betas**2
print(ombc)

print("------------------fac=qV0*omeag/mc2/c/(bg)**3")
fac = qV0*omega/(m0c2*c*(gammas*betas)**3)
print(fac)
print(Ws*betas**(-0.5)*fac*ombc/gammas**4/m0c2)
print()

print("------------------------T'(k)")
k, Dz, a, b = symbols('k Dz a, b')
Tpks = -(2*sin(k*Dz)/(k*(2*Dz+2/3*b*Dz**3)))*((1+3*b*Dz**2-6*b/k**2)/k-Dz*cot(k*Dz)*(1+b*Dz**2-6*b/k**2))
print(Tpks)

print("-----------------------T''(k)")
Tppks = diff(Tpks,k)
print (simplify(Tppks))

print("-----------------------S'(k)")
Spks = 2*a*sin(k*Dz)/(k*(2*Dz+2/3*b*Dz**2))*(Dz**2-2/k**2+Dz*cot(k*Dz)*2/k)
print(Spks)

print("---------------------------S''(k)")
Sppks = diff(Spks,k)
print(simplify(Sppks))

