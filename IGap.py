from abc import ABC,abstractmethod

class IGap(ABC):
    @abstractmethod
    def configure(self,**kwargs): pass
    @abstractmethod
    def map(self,i_track): pass
    @abstractmethod
    def toString(self): pass
    @abstractmethod
    def isAccelerating(self): pass
    @abstractmethod
    def adjust_energy(self, tkin): pass
    @abstractmethod
    def waccept(self,**kwargs): pass
    @abstractmethod
    def register(self): pass
    @abstractmethod
    def accept(self): pass
