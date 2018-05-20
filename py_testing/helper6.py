import inspect

class Aclass(object):
    def __init__(self,magic='?'):
        self.name  = "Aclass"
        self.label = "Aclass"
        self.magic = magic
        

class Bclass(Aclass):
    """ Documentation of Bclass """
    def __init__(self):
        super().__init__('abacadraba')
        # super().__init__()
        self.name = "Bclass"


b=Bclass()

print('{}\n'.format(b))
print(''.join('{}\n'.format(el) for el in inspect.getmembers(b)))
print('MRO:\n'+''.join('{}\n'.format(el) for el in Bclass.__mro__))
