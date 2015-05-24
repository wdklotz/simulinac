#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
# from setup import Phys,Beam,Proton,objprnt
# from math import pi,sqrt,sin
import yaml

def test0():
    wfl= []
    file=open('template.yml','r')
    wfl= yaml.load(file)
    print(yaml.dump(wfl,default_flow_style=True))
    for i,v in iter(wfl.items()):
        print(i,' =\t',v)
    seg = wfl['segment']
    print(seg)
    print('=== segment ===')
    for i in seg:
        print(i)
    lattice = wfl['lattice']
    print('=== lattice ===')
    for l in lattice:
        for i in l:
            print(i)
def test1():
    fin= open('template.yml','r')
    data_in= yaml.load(fin)

    flags= data_in['flags']
    print('\nflags=\t',flags)

    parameters= data_in['parameters']
    print('\nparameters=\t',parameters)

    elements= data_in['elements']
    print('\nelements=\t',elements)

    segments= data_in['segments']
    print('\nsegments=\t',segments)

    lattice= data_in['lattice']
    print('\nlattice=\t',lattice)
#--------*--------*--------*--------*--------*--------*--------*--------*--------*
if __name__ == '__main__':
#     test0()
    test1()
