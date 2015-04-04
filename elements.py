# -*- coding: utf-8 -*-
import setup
import numpy as np
from math import sqrt, sinh, cosh, sin, cos, fabs, tan, floor, modf

class Matrix(object):
    _dim = 4
    def __init__(self):
        self.matrix=np.eye(Matrix._dim)    ## 4x4 unit matrix
        self.label='I'
        self.length=0.         ## default zero length!
        self.slice_min = 0.01  ## minimal slice length
        self.viseo = 0.
    #-----------------------
    def out(self):
        print(self.label)
        print(self.matrix)
    #-----------------------
    def __mul__(self,other):
        product=np.einsum('ij,jk',self.matrix,other.matrix)
        res=Matrix()
        res.label=self.label+'*'+other.label
        res.matrix=product
        return res
    #-----------------------
    def reverse(self):
        res=Matrix()
        for i in range(4):
            for k in range(Matrix._dim):
                res.matrix[i][k] = self.matrix[i][k]
        res.matrix[0][0] = self.matrix[1][1]
        res.matrix[1][1] = self.matrix[0][0]
        res.matrix[2][2] = self.matrix[3][3]
        res.matrix[3][3] = self.matrix[2][2]
        res.label = '('+self.label+')r'
        return res
    #-----------------------
    def trace(self):
        return self.tracex()+self.tracey
    def tracex(self):
        res = 0.
        for i in range(2):
            res += self.matrix[i,i]
        return res
    def tracey(self):
        res = 0.
        for i in range(2,4):
            res += self.matrix[i,i]
        return res
    #-----------------------
    def shorten(self,length=0.):
        return Matrix()
    #-----------------------
    def step_through(self,anz=10):
        """
        Step through an element(The central nontrivial function).
        Default is 10 steps/element.
        Minimal step size is self.slice_min.
        """
        if self.length == 0.0:   ## zero length element (like wedge)
            yield self
        else:                    ## adjust step size and nboff steps
            step = self.length/anz
            if step < self.slice_min:
                step = self.slice_min
#            fanz = self.length/step
#            anz  = floor(fanz)
#            rest = self.length - anz*step
            (rest,fanz) = modf(self.length/step)
            anz = int(fanz)
            rest = self.length * rest
#            print('fanz anz, step, rest',fanz,anz,step,rest)
            mx = self.shorten(step)
            if fabs(rest) > 1.e-12:
                mr = self.shorten(rest)
            for i in range(anz+1):
                if i == anz:
                    if fabs(rest) < 1.e-12:
                        break
                    else:
                        yield mr
                else:
                    yield mx
    #-----------------------
    def BetaMatrix(self):
        m11 =self.matrix[0,0];  m12 =self.matrix[0,1]
        m21 =self.matrix[1,0];  m22 =self.matrix[1,1]
        n11 =self.matrix[2,2];  n12 =self.matrix[2,3]
        n21 =self.matrix[3,2];  n22 =self.matrix[3,3]
        m_beta = np.array([
            [ m11*m11, -2.*m11*m12,           m12*m12,   0., 0., 0.],
            [-m11*m21,     m11*m22+m12*m21,  -m22*m12,   0., 0., 0.],
            [ m21*m21, -2.*m22*m21,           m22*m22,   0., 0., 0.],
            [ 0., 0., 0., n11*n11, -2.*n11*n12,           n12*n12],
            [ 0., 0., 0.,-n11*n21,     n11*n22+n12*n21,  -n22*n12],
            [ 0., 0., 0., n21*n21, -2.*n22*n21,           n22*n22]
            ])
        return m_beta

class I(Matrix):           ## unity Matrix (an alias to Matrix class)
    def __init__(self):
        super(I,self).__init__()
        
class Test(Matrix):
    def __init__(self,a,b,c,d,label='test'):
        self.matrix=np.array([[a,b,0.,0.],
                              [c,d,0.,0.],
                              [0.,0.,1.,2.],
                              [0.,0.,-1.,1.]])
        self.label=label

class D(Matrix):     ## drift space
    def __init__(self,length=0.,label='D'):    
        super(D,self).__init__()
        self.label=label
        self.length=length     ## hard edge length [m]
        self.matrix[0][1]=self.matrix[2][3]=self.length
    def shorten(self,l=0.):
        return D(length=l,label=self.label)
        

class QF(D):    ## focusing quad
    def __init__(self,k0=0.,length=0.,label='QF'):    
        super(QF,self).__init__(length=length,label=label)
        self.k0=k0         ## energy independant Quad strength [m**-2]
        self.matrix=self._mx()
        self.viseo = +0.5
    def shorten(self,l=0.):
        return QF(k0=self.k0,length=l,label=self.label)
    def _mx(self):
        kwurz=sqrt(self.k0)
        phi=self.length*kwurz
        ## focusing
        cf  =cos(phi)
        sf  =sin(phi)/kwurz
        cfp =-kwurz*sin(phi)
        sfp =cf
        ## defocusing
        cd  =cosh(phi)
        sd  =sinh(phi)/kwurz
        cdp =kwurz*sinh(phi)
        sdp =cd
        ## 4x4 matrix
        if (isinstance(self,QF) and (isinstance(self,QD)==False)):
            mq=np.array([[cf,sf,0.,0.],[cfp,sfp,0.,0.],[0.,0.,cd,sd],[0.,0.,cdp,sdp]])  ## QF
        elif isinstance(self,QD) :
            mq=np.array([[cd,sd,0.,0.],[cdp,sdp,0.,0.],[0.,0.,cf,sf],[0.,0.,cfp,sfp]])  ## QD
        else:
            raise RuntimeError('QF._mx: neither QF nor QD! should never happen!')
        return mq

class QD(QF):                 ## defocusing quad
    def __init__(self,k0=0.,length=0.,label='QD'):
        super(QD,self).__init__(k0,length,label)
        self.viseo = -0.5
    def shorten(self,l=0.):
        return QD(k0=self.k0,length=l,label=self.label)

class SD(D):   ## sector bending dipole in x-plane
    def __init__(self,radius=0.,length=0.,label='SB'):
        super(SD,self).__init__(length=length,label=label)
        self.radius = radius
        self.matrix=self._mx()
        self.viseo = 0.25
    def shorten(self,l=0.):
        return SD(radius=self.radius,length=l,label=self.label)
    def _mx(self):
        kwurz=1./self.radius
        phi=self.length*kwurz
        ## x-plane
        cf  = cos(phi)
        sf  = sin(phi)/kwurz
        cfp =-kwurz*sin(phi)
        sfp = cf
        ## y-plane
        cd  = 1.
        sd  = phi/kwurz
        cdp = 0.0
        sdp = cd
        ## 4x4 matrix
        ms=np.array([[cf,sf,0.,0.],[cfp,sfp,0.,0.],[0.,0.,cd,sd],[0.,0.,cdp,sdp]])   ## sector
        return ms

class WD(D):   ## wedge of rectangular bending dipole in x-plane
    def __init__(self,length=0.,radius=0.,label='WD'):
        super(WD,self).__init__(label=label)
        self.radius = radius
        kwurz=1./self.radius
        psi=0.5*length*kwurz        ## Kantenwinkel
        ckp=kwurz*tan(psi)
        ## 4x4 matrix
        mw=np.array([[1.,0.,0.,0.],[ckp,1.,0.,0.], [0.,0.,1.,0.],[0.,0.,-ckp,1.]])   ## wedge
        self.matrix=mw
#####################################################################
def k0(gradient=0.,beta=0.,energy=0.):   ## helper function for tests
    """
        quad strength as function of energy and gradient
        gradient in [Tesla/m]
        energy in [Gev]
        beta [v/c]
    """
    if (gradient != 0. and energy != 0. and beta !=0.):
        return 0.2998*gradient/(beta*energy)             
    else:
        raise RuntimeError('zero gradient or energy or beta in quad strength!')

def test0():
    print('trivial test 0 ...')
    a=Test(1,2,3,4,'a')
    a.out()
    b=Test(1,0,3,4,'b')
    b.out()
    (a*b).out()
    (b*a).out()
    
def test1():
    print('trivial test 1 ...')
    i1=Matrix()
    i2=i1*i1
    i1.out()
    i2.out()
    
def test2():
    print('trivial test 2 ...')
    i1=Matrix()
    d1=D(10.,'D1')
    d1.out()
    (d1*d1).out()
    (d1*d1*d1).out()
    (d1*i1*d1).out()
    d2=D(90,'D2')
    d2.out()
    (d1*d2).out()

def test3():
    print('test product of Matrix class ...')
    gradient =1.
    beta     =0.5
    energy   =0.2
    print('gradient[Tesla/m] {:.3f}; beta[v/c] {:.3f}; energy[Gev] {:.3f}'.format(gradient,beta,energy))
    k=k0(gradient=gradient,energy=energy,beta=beta)
    qf=QF(k0=k,length=1.)
    qf.out()
    ## test product of Matrix class
    qd=QD(k0=k,length=1.)
    qd.out()
    (qf*qd).out()

def test4():
    print('test shortening of elements ...')
    gradient =1.
    beta     =0.5
    energy   =0.2
    k=k0(gradient=gradient,energy=energy,beta=beta)
    ## elements
    d10=D(10.,'d10')
    d10.out()
    (d10.shorten(1.e-2)).out()
    
    qf=QF(k0=k,length=1.)
    qf05=qf.shorten(0.2)
    (qf05*qf05*qf05*qf05*qf05).out()
    qf.out()
    
    qd=QD(k0=k,length=1.)
    qd05=qd.shorten(0.2)
    (qd05*qd05*qd05*qd05*qd05).out()
    qd.out()
    
    sd=SD(radius=10.,length=2.)
    sd05=sd.shorten(0.4)
    (sd05*sd05*sd05*sd05*sd05).out()
    sd.out()

def test5():
    print("K.Wille's Beispiel auf pp. 112-113")
    kqf=  setup.ex_wille()['k_quad_f']
    lqf=  setup.ex_wille()['length_quad_f']
    kqd=  setup.ex_wille()['k_quad_d']
    lqd=  setup.ex_wille()['length_quad_d']
    rhob= setup.ex_wille()['beding_radius']
    lb=   setup.ex_wille()['dipole_length']
    ld=   setup.ex_wille()['drift_length']
    ## elements
    mqf=QF(kqf,lqf,'QF')
    mqd=QD(kqd,lqd,'QD')
    mb=SD(rhob,lb,'B')
    mw=WD(lb,rhob)
    md=D(ld)
    ## test matrix multiplication
    mz=I()
    mz=mz *mqf
    mz=mz *md
    mz=mz *mw
    mz=mz *mb
    mz=mz *mw
    mz=mz *md
    mz=mz *mqd
    mz=mz *md
    mz=mz *mw
    mz=mz *mb
    mz=mz *mw
    mz=mz *md
    mz=mz *mqf
    mz.out()

def test6():
    print('test step_through elements ...')
    kqf=  setup.ex_wille()['k_quad_f']
    lqf=  setup.ex_wille()['length_quad_f']
    kqd=  setup.ex_wille()['k_quad_d']
    lqd=  setup.ex_wille()['length_quad_d']
    rhob= setup.ex_wille()['beding_radius']
    lb=   setup.ex_wille()['dipole_length']
    ld=   setup.ex_wille()['drift_length']

    ## elements
    mqf=QF(kqf,lqf,'QF')
    mqd=QD(kqd,lqd,'QD')
    mb=SD(rhob,lb,'B')
    mw=WD(lb,rhob)
    md=D(ld)
    
    ## test step_through elements ...
    list=[mqf,mqd,mb,mw,md]
    # list=[mqf]
    for m_anfang in list:
        m_end=I()
        print('======================================')
        for count,mi in enumerate(m_anfang.step_through()):
            print('step ',count+1,end='  ')
            mi.out()
            m_end=m_end*mi
        m_end.out()
        m_anfang.out()

def test7():
    print('======================================')
    print('test Rechteckmagnet...')
    rhob= setup.ex_wille()['beding_radius'] 
    lb=   setup.ex_wille()['dipole_length']
    mb=SD(radius=rhob,length=lb,label='B')
    mw=WD(length=lb,radius=rhob,label='W')
    mr=mw*mb*mw
    mw.out()
    mb.out()
    mr.out()
##################################################################
if __name__ == '__main__':
    test0()
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
    test7()
