class Base:
    def __init__(self):
        self.name = self.__class__.__name__
        self.var1 = None
        self.var2 = None
    # def pn(self):
        # print(self.name)
    def pv(self):
        print(self.name,self.var1, self.var2)

class Subclass(Base):
    def __init__(self,var2):
        Base.__init__(self)
        self.var2 = var2

b=Base()
s=Subclass("var2i")
print(b.__dict__)
print(s.__dict__)
print("======A")
b.pv()
s.pv()
print("======B")
b.var1 = "b_var1"   # set a Base property
s.var2 = "s_var2"   # set a Subclass property
b.pv()
s.pv()
b.pv()
print("======C")
s.var1 = "s_var1"   # set a inherited base property
b.pv()
s.pv()
print("======D")
print(b.__dict__)
print(s.__dict__)
