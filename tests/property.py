class Base(object):
    def __init__(self):
        self._number=999
 
@property
def name(self):
    '''I am the name property'''
    return self._name
@name.setter
def name(self,value):
    self._name=value

def shownumber(self):
    print('my number: {}'.format(self._number))

Base.name   = name
Base.number = shownumber
help(Base.name)

base = Base()
base.name = 'xxxxx'
print(base.name)
print(base.__dict__)
base.number()

