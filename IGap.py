from abc import ABC,abstractmethod

class IGap(ABC):
    @abstractmethod
    def configure(self,**kwargs): pass

    @abstractmethod
    def values_at_exit(self): pass

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
    def values_at_exit(self): pass

    @abstractmethod
    def register_mapper(self,obj): pass

    @abstractmethod
    def accept_register(self): pass
