import os
import numpy as np

template = 'WorkOnArc3.yml'
simu_input = 'WOA2.yml'

# seq = np.arange(1.37735, 1.37745, 0.00001)
# seq  = np.arange(1.15, 1.25, 0.01)
# seq1 = np.arange(1.35, 1.45, 0.01)

# for v in seq:
#     for v1 in seq1:
#         command = f'm4 -D _Gf={v:.6f} -D _Gd={v1:.6f} {template} > {simu_input}'
#         print(command)
#         os.system(command)
#         command = f'python simu.py --file {simu_input}'
#         os.system(command)
#         print(f'==================>_Gf = {v:.5f} _Gd = {v1:.5f}')

seq  = np.arange(6.37, 2*6.37, 0.03)

for v in seq:
    command = f'm4 -D _RHO={v:.6f} {template} > {simu_input}'
    print(command)
    os.system(command)
    command = f'python simu.py --file {simu_input}'
    os.system(command)
    print(f'==================>_RHO = {v:.5f}' )
