import sys
from math import sqrt,fabs,acos,degrees
from numpy import linalg as LA
import numpy as NP
from copy import copy
import warnings

from lattice import Lattice
# from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO
from setutil import wille,PARAMS,FLAGS,SUMMARY,printv,DEB,sigmas, objprnt, Ktw, Ktp
# from setutil import Twiss, Functions
# import elements as ELM
# import TTFG as TTF
# from sigma import Sigma

DEBUG_MODULE = DEB.get('OFF')

class Section(Lattice):
    """
    A Section is part of a Lattice
    """
    def __init__(self) -> None:
        super().__init__()

