import sys
sys.path.insert(0,'..')

from DynacG import Picht as Picht1
import numpy as NP

T=50.09747033226243
# T=70
gamma = 1.+T/938.272

r=[1,1]
rp=[0,0]
DgDz=0.02
R,Rp=Picht1(gamma,r,rp,DgDz)
print(r,rp,'\n',R,Rp)
r,rp=Picht1(gamma,R,Rp,DgDz,back=True)
print(R,Rp,'\n',r,rp)
print()
r=[0,0]
rp=[1,1]
R,Rp=Picht1(gamma,r,rp,DgDz)
print(r,rp,'\n',R,Rp)
r,rp=Picht1(gamma,R,Rp,DgDz,back=True)
print(R,Rp,'\n',r,rp)
print()
r=[1,1]
rp=[1,1]
R,Rp=Picht1(gamma,r,rp,DgDz)
print(r,rp,'\n',R,Rp)
r,rp=Picht1(gamma,R,Rp,DgDz,back=True)
print(R,Rp,'\n',r,rp)
print()
r=[1,0]
rp=[0,1]
R,Rp=Picht1(gamma,r,rp,DgDz)
print(r,rp,'\n',R,Rp)
r,rp=Picht1(gamma,R,Rp,DgDz,back=True)
print(R,Rp,'\n',r,rp)
print()
r=[0,1]
rp=[1,0]
R,Rp=Picht1(gamma,r,rp,DgDz)
print(r,rp,'\n',R,Rp)
r,rp=Picht1(gamma,R,Rp,DgDz,back=True)
print(R,Rp,'\n',r,rp)
