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
