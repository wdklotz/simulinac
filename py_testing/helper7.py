class A():
    """testing and understanding __init_subclass__"""
    def __init_subclass__(cls,**kwargs):
        print(cls,kwargs)
        for k in kwargs:
            setattr(cls,k,kwargs[k])
        def out(cls):
            print('p1={}, p2={}, p3={}'.format(cls.p1,cls.p2,cls.p3))
        setattr(cls,"out",out)


class B(A,p1=1,p2=2,p3=3):
    """class B(A,**kwargs)"""
    pass

print(A.__dict__)
print(B.__dict__)   
b=B()
print('p1={}, p2={}, p3={}'.format(b.p1,b.p2,b.p3))
b.p1=11;b.p2=22;b.p3=33
print('p1={}, p2={}, p3={}'.format(b.p1,b.p2,b.p3))
b.p1=111;b.p2=222;b.p3=333
b.out()

class Philosopher:
    def __init_subclass__(cls, default_name, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.default_name = default_name

class AustralianPhilosopher(Philosopher, default_name="Bruce"):
    pass

bruce = AustralianPhilosopher()
print(bruce.default_name)
