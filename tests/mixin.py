class Base(object):
    def __init__(self):
        self.number=999
    def attributes(self):
        print('Base.number= ',self.number)

class Mixin1(object):
    def set(self,name):
        self.name=name.upper()
    def get(self):
        return self.name

class Mixin2(object):
    def set(self,name):
        self.name='{}:{}'.format(name,name)
    def get(self):
        return self.name

class Mixin3(object):
    def set(self,name):
        self.name=('{}:{}'.format(name,name)).upper()
    def get(self):
        return self.name

class NamedBase(Mixin1,Mixin2,Base):
    def __init__(self):
        super().__init__()
    def who(self):
        print('I am {} '.format(self.name))

base = NamedBase()
base.set('NamedBase')
print('bases: ',base.__class__.__bases__)
print('subclasses: ',base.__class__.__subclasses__())
print('dict: ',base.__dict__)
base.who()
base.attributes()
assert isinstance(base,Base)
assert isinstance(base,NamedBase)
assert isinstance(base,Mixin1)
assert isinstance(base,Mixin2)
assert issubclass(NamedBase,Mixin1)
# assert isinstance(base,Mixin3)

import inspect
for i, cls in enumerate(inspect.getmro(NamedBase)):
    print(i,cls)