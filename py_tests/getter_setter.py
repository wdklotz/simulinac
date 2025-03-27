class Common():
    def __init__(self):
        self.a=None
        self.b=None
        self.c=None
        self.d=None
    
class Top():
    def __init__(self,common,sub=None):
        self.com=common
        self.sub=sub

    # readables
    @property
    def a(self): return self.com.a
    @property
    def b(self): return self.com.b
    @property
    def c(self): return self.com.c
    @property
    def d(self): return self.com.d

    # writables
    @c.setter
    def c(self,v): self.com.c=v
    @d.setter
    def d(self,v): self.com.d=v

class Sub():
    def __init__(self,common,top):
        self.com=common
        self.top=top

    # readables
    @property
    def a(self): return self.com.a
    @property
    def b(self): return self.com.b
    @property
    def c(self): return self.com.c
    @property
    def d(self): return self.com.d

    # writables
    @c.setter
    def c(self,v): self.com.c=v
    @d.setter
    def d(self,v): self.com.d=v

common=Common()
top=Top(common)
sub=Sub(common,top)
top.sub=sub

print(top.a)   
print(sub.a)

sub.c=99
print(top.c)

top.d=4
print(sub.d)

top.com.a=99
print(sub.a)

top.a=99
