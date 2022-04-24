import pprint
pp =pprint.PrettyPrinter()

class Base(object):
    def __init__(self):
        self._number=999
 
@property
def name(self):
    '''I am the name property'''
    return self._name
@name.setter
def name(self,value):
    self._name = value

def shownumber(self):
    print('method added to class prints my number: {}'.format(self._number))

print('-------- adding a member to a class --------')
print('before: ')
pp.pprint(Base.__dict__)
Base.shownumber = shownumber
print('after: ')
pp.pprint(Base.__dict__)

print('-------- adding an attribute to a class --------')
print('before: ')
pp.pprint(Base.__dict__)
Base.attribute = 'ATTRIBUTE'
print('after: ')
pp.pprint(Base.__dict__)

print('-------- adding a property to a class --------')
pp.pprint(name)
print('before: ')
pp.pprint(Base.__dict__)
Base.name   = name
print('after: ')
pp.pprint(Base.__dict__)
help(Base.name)

base = Base()
print('-------- setting the object property --------')
print('before: ')
pp.pprint(base.__dict__)
base.name = 'NAME'
print('after: ')
pp.pprint(base.__dict__)
print('-------- getting the object property --------')
print(base.name)

