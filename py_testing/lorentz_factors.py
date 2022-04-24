import math as M
p_mass = 0.9382720813 # GeV
tk = 5000.e-3            # GeV (5 MeV)

E = p_mass + tk
# pc = M.sqrt(E**2-p_mass**2)
# print('pc [GeV/c] {}'.format(pc))
# brho = 3.3356*pc
# print('brho [T*m] {}'.format(brho))

gamma = 1.+tk/p_mass
gb = M.sqrt(gamma**2-1)
beta = gb/gamma
pc=beta*E
print('pc=beta*E   [GeV/c] {}'.format(pc))
pc = gb*p_mass
print('pc=g*b*m0c2 [Gev/c] {}'.format(pc))
brho = 3.1297 * gb
print('brho        [T*m]   {}'.format(brho))
print('gamma {}\nbeta  {}\ng*b   {}'.format(gamma,beta,gb))
print()

b2 = 1.
print('b2 {} b2(scaled) {}'.format(b2,b2/brho))
