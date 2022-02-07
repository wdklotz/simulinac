class A:
    def __init__(self):
        self.values = {}
    def set(self,k,v):
        self.values[k] = v
    def __setitem__(self,k,v):
        self.values[k] = v
    def __getitem__(self,k):
        return self.values[k]

aclass = A()
kv=('a',1)
aclass.set(*kv)        # use set(self,...)
aclass['b'] = 2        # use __setitem__(self,...)
print(aclass.__dict__) # use builtin __dict__
print(aclass['b'])     # use __getitem__(self,...)
 