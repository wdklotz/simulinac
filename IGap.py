#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
__version__='v11.0.3'
"""
Copyright 2015 Wolf-Dieter Klotz <wdklotz@gmail.com>
This file is part of the SIMULINAC code

    SIMULINAC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    SIMULINAC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SIMULINAC.  If not, see <http://www.gnu.org/licenses/>.

"""

from abc import ABC,abstractmethod

class IGap(ABC):
    @abstractmethod
    def accelerating(self): pass

    @abstractmethod
    def adjust_energy(self, tkin): pass

    @abstractmethod
    def configure(self,**kwargs): pass

    @abstractmethod
    def map(self,i_track): pass

    @abstractmethod
    def register(self,obj): pass

    @abstractmethod
    def toString(self): pass

    @abstractmethod
    def waccept(self,**kwargs): pass
