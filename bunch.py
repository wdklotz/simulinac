#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
"""
Copyright 2015 Wolf-Dieter Klotz <wdklotz@gmail.com>
This file is part of the SIMULINAC code

    SIMULINAC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    SIMULINAC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SIMULINAC.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as NP
from math import sqrt
import matplotlib.pyplot as plt

from setutil import DEBUG, Proton, tblprnt, K, sigmas, PARAMS
from trackPlot import histPlot, poincarePlot
from Dictionary import DictObject

# DEBUG
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass

class Tpoint(object):
    """ 
        A track point is an NP.array of 10 coordinates, 
        i.e. (0=x, 1=x', 2=y, 3=y', 4=z, 5=z', 6=T, 1, 8=s, 1) 
    """
    def __init__(self, point = NP.array([0,0,0,0,0,0,0,1,0,1])):
        self.point = point
    def __call__(self):
        return self.point
    def as_str(self):
        s = 'x={:10.03e} x\'={:10.03e} y={:10.03e} y\'={:10.03e} z={:10.03e} z\'={:10.03e} T={:7.02f}  S={:7.02f}'.format(self.point[K.x],self.point[K.xp],self.point[K.y],self.point[K.yp],self.point[K.z],self.point[K.zp],self.point[K.T],self.point[K.S])
        return s
    

class Track(DictObject,object):
    """
        A Track is a dictionary
        A Track is a container (tuple) of positions (Tpoint). 
    """
    def __init__(self):
        DictObject.__init__(self)
        self.points = []
    def __getitem__(self,n):        # iterator over points in Track
        return self.points[n]

    def getpoints(self):            # points in Track
        return self.points
    def nbofpoints(self):           # nbof points in Track
        return len(self.points)
    def addpoint(self,point):       # add a point to Track
        if self.nbofpoints() == 0:
            last = 0
        else:
            last, p = self.points[-1]
        self.points.append((last+1,point))
    def as_table(self):
        tblheadr = ['    x',"    x'",'    y',"    y'",'    z',"    z'",'  tkin','    s']
        # tblheadr = ['x',"x'",'y',"y'",'z',"z'",'tkin','s']
        tblrows =[]
        for n,point in iter(self):
            tblrow = [
                '{:8.3f}'.format(point()[K.x]),
                '{:8.3f}'.format(point()[K.xp]),
                '{:8.3f}'.format(point()[K.y]),
                '{:8.3f}'.format(point()[K.yp]),
                '{:8.3f}'.format(point()[K.z]),
                '{:8.3f}'.format(point()[K.zp]),
                '{:8.3f}'.format(point()[K.T]),
                '{:8.3f}'.format(point()[K.S]),
                ]
            tblrows.append(tblrow)
        return tblprnt(tblheadr,tblrows)
    def as_str(self):
        str = ''
        for n,p in iter(self):
            str += p.as_str()+'\n'
        return str

class Bunch(DictObject,object):
    """
        A Bunch is a dictionary.
        A Bunch is a container (tuple) of particles
    """
    def __init__(self):
        DictObject.__init__(self)
        self.particles = []
    def __getitem__(self,n):            # iterator over particles in bunch
        return self.particles[n]

    def getparticles(self):             # particles in bunch
        return self.particles
    def nbofparticles(self):            # nbof particles in bunch
        return len(self.particles)
    def addparticle(self,particle):     # add particle to bunch
        if self.nbofparticles() == 0:
            last = 0
        else:
            last, p = self.particles[-1]
        self.particles.append((last+1,particle))
    # @property
    # def distfactory(self):  # distribution factory
    #     return self['disttype']
    # @distfactory.setter
    # def distfactory(self,value):
    #     self['disttype'] = value

   ##   def populate_phase_space(self):
    #     """ create the distribution """
    #     self.tracklist = self.disttype(self)

class EmitContour(object):
    """
        Generates a bunch with particles of same emittance
    """
    def emittanceContourPoint(self, x, alfa, beta, emit):
        gamma = (1.+ alfa**2)/beta
        a = beta
        b = 2.*alfa*x
        c = gamma*x**2-emit
        d = sqrt(b**2-4.*a*c)
        y = (-b+d)/(2.*a)   
        p1 = (x,y)     # upper half-plane
        p2 = (-x,-y)   # upper half-plane 
        return (p1,p2)

    def __init__(self,nTracks,random=False):
        sigx,sigxp = sigmas(PARAMS['alfax_i'],PARAMS['betax_i'],PARAMS['emitx_i'])
        if random:
            Xrand = sigx*(2.*NP.random.random_sample((nTracks,))-1.)
        else:
            Xrand = NP.linspace(-sigx*(1.-1.e-3),sigx*(1.-1.e-3),nTracks)
        self.X=[]; self.XP=[]
        for x in Xrand:
            points = self.emittanceContourPoint(x,PARAMS['alfax_i'],PARAMS['betax_i'],PARAMS['emitx_i'])
            self.X.append( points[0][0])
            self.XP.append(points[0][1])
            self.X.append( points[1][0])
            self.XP.append(points[1][1])
        sigy,sigyp = sigmas(PARAMS['alfay_i'],PARAMS['betay_i'],PARAMS['emity_i'])
        if random:
            Yrand = sigy*(2.*NP.random.random_sample((nTracks,))-1.)
        else:
            Yrand = NP.linspace(-sigy+1.e-5,sigy-1.e-5,nTracks)
        self.Y=[]; self.YP=[]
        for y in Yrand:
            points = self.emittanceContourPoint(y,PARAMS['alfay_i'],PARAMS['betay_i'],PARAMS['emity_i'])
            self.Y.append( points[0][0])
            self.YP.append(points[0][1])
            self.Y.append( points[1][0])
            self.YP.append(points[1][1])
        tkin = PARAMS['sollteilchen'].tkin  #energy at entrance
        track = Track()
        for i in range(2*nTracks):
            start       = NP.array([ 0., 0., 0., 0., 0., 0., tkin, 1., 0., 1.])
            start[K.x]  = self.X[i]
            start[K.xp] = self.XP[i]
            start[K.y]  = self.Y[i]
            start[K.yp] = self.YP[i]
            point = Tpoint(start)
            track.addpoint(point)
        self.track = track

def Gauss1D(params):
    """ generates a bunch with 1D gaussian distribution """
    sigx       = params['sigx']
    sigxp      = params['sigxp']
    sigy       = params['sigy']
    sigyp      = params['sigyp']
    sigz       = params['sigz']
    sigzp      = params['sigzp']
    nbtracks   = params['nbtracks']
    coord_mask = params['coord_mask']
    tkin       = params['tkin']

    X          = sigx  * NP.random.randn(nbtracks)    # gauss distribution X
    XP         = sigxp * NP.random.randn(nbtracks)
    Y          = sigy  * NP.random.randn(nbtracks)
    YP         = sigyp * NP.random.randn(nbtracks)
    Z          = sigz  * NP.random.randn(nbtracks)
    ZP         = sigzp * NP.random.randn(nbtracks)

    tracklist=[]        # all tracks in a bunch
    for i in range(nbtracks):
        start=NP.array([ 0., 0., 0., 0., 0., 0., tkin, 1., 0., 1.])
        # initial setting for each coordinate
        if coord_mask[K.x]:
            start[K.x]  = X[i]
        if coord_mask[K.xp]:
            start[K.xp] = XP[i]
        if coord_mask[K.y]:
            start[K.y]  = Y[i]
        if coord_mask[K.yp]:
            start[K.yp] = YP[i]
        if coord_mask[K.z]:
            start[K.z]  = Z[i]
        if coord_mask[K.zp]:
            start[K.zp] = ZP[i]
        track = Track()
        tp = Tpoint(point=start)
        track.addpoint(tp)
        tracklist.append(track)
    return tracklist

def test0():
    print('-----------------------------------------Test0---')
    bunch = Bunch()
    print('nbofparticles: ', bunch.nbofparticles())
    # populate bunch
    for i in range(3):
        p = Proton()
        bunch.addparticle(p)
    print('nbofparticles: ', bunch.nbofparticles())
    allparticles = bunch.getparticles()
    print('particles: ', allparticles)
    print('bunch.__dict__: ', bunch.__dict__)
    # loop particles in bunch
    for particle in iter(bunch):
        print(particle)
    # populate track
    track = Track()
    for i in range(3):
        point = Tpoint(NP.array([i,0,i,0,i,0,0,1,0,1]))
        track.addpoint(point)
    print('nbofpoints: ',track.nbofpoints())
    allpoints = track.getpoints()
    print('points: ',allpoints)
    print('track.__dict__: ',track.__dict__)
    # loop points in track
    for n,point in iter(track):
        print(n,point())
    # link track with particle
    for n,particle in iter(bunch):
        particle['track'] = track
    # show
    particles = bunch.getparticles()
    n, last = particles[-1]
    print('last particle track: ', last['track'])
    print("last['track'].__dict__; ",last['track'].__dict__)
    
    print(last['track'].as_table())

def test1():
    # example data
    mu    = 0     # mean of distribution
    sigma = 1     # standard deviation of distribution
    print('-----------------------------------------Test1---')
    # NP.random.randn returns a sample (or samples) from the “standard normal” distribution.
    x   = mu + sigma * NP.random.randn(4000)
    fig = plt.figure('test1: figure')
    histPlot(x,mu,sigma)
    figures.append(fig)

def test2():
    print('-----------------------------------------Test2---')
    N = 20000
    alfax = 1
    betax = 1.
    emitx = 1.e-3
    sigma, sigmap = sigmas(alfax, betax, emitx)
    x      = sigma  * NP.random.randn(N)
    xp     = sigmap * NP.random.randn(N)
    DEBUG_OFF('x: ', x)
    DEBUG_OFF('x\':',xp)

    fig1 = plt.figure('test2:figure 1')
    h1   = plt.subplot2grid((2,1),(0,0))
    histPlot(x,0.,sigma)
    h2   = plt.subplot2grid((2,1),(1,0))
    histPlot(xp,0.,sigmap)
    figures.append(fig1)

    fig2 = plt.figure('test2:figure 2')
    good = (x, xp)
    poincarePlot(good, 'x-x\'', (0.1,0.1), projections=(1,1))
    figures.append(fig2)

def test3(filepath):
    print('-----------------------------------------Test3---')
    from lattice_generator import parse_and_fabric
    N = 200
    lattice = parse_and_fabric(filepath)
    contour = EmitContour(N, random=True)
    track   = contour.track
    X  = [x()[K.x] for n,x in iter(track)]
    XP = [x()[K.xp] for n,x in iter(track)]
    Y  = [x()[K.y] for n,x in iter(track)]
    YP = [x()[K.yp] for n,x in iter(track)]
    fig = plt.figure('test2:figure')
    plt.scatter(X,XP,s=0.1)
    plt.scatter(Y,YP,color='red',s=0.1)
    figures.append(fig)

def test4():
    print('-----------------------------------------Test4---')
    params = dict(
    sigx       = 6,
    sigxp      = 1,
    sigy       = 1,
    sigyp      = 6,
    sigz       = 3,
    sigzp      = 10,
    nbtracks   = 30000,
    coord_mask = NP.array([1,1,1,1,1,1,0,1,0,1]),
    tkin       = 70
    )
    
    tracklist = Gauss1D(params)
    xaxis = []
    yaxis = []
    for track in tracklist:
        for n,point in iter(track):
            X  = point()[K.x]
            XP = point()[K.xp]
            xaxis.append(X)
            yaxis.append(XP)
    fig = plt.figure("test4:figure")
    plt.scatter(xaxis,yaxis,s=0.1)
    figures.append(fig)
            
if __name__ == '__main__':
    figures = []
    test0()
    test1()
    test2()
    test3('yml/work.yml')
    test4()
    [plt.show(fig) for fig in figures]
