#!/Users/klotz/SIMULINAC_env/bin/python
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
import numpy as np
from copy import copy

from setutil import PARAMS,DEBUG
from elements import MDIM,XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
import elements as ELM

## DEBUG MODULE
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_MODULE = DEBUG_OFF

## Track class
class Track(object):    #is an ordered list of track-points. A track-point is an array of MDIM coordinates.

    def __init__(self, particle_number=0, start=None):
        self.track_points = start
        self.particle_number = particle_number
        self.nb_points_per_track = 1

    def nb_points(self):
        return self.nb_points_per_track

    def append(self,new):
        self.track_points = np.append(self.track_points,new)
        self.nb_points_per_track +=1

    def points(self):
        return self.track_points.reshape(self.nb_points_per_track,MDIM)

    def point_at(self,n):
        return self.points()[n]

    def first(self):
        first = self.point_at(0)
        return first

    def last(self):
        last = self.point_at(-1)
        return last

    def first_str(self):
        return Track.string(self.first())

    def last_str(self):
        return Track.string(self.last())

    def points_str(self):
        points = self.points()
        str = ''
        for p in points:
            str += Track.string(p)
        return str

    def string(p):   #single point to string
        s = 'x={:.3e} x\'={:.3e} y={:.3e} y\'={:.3e} z={:.3e} z\'={:.3e}  tk={:.5f} s={:.3f} '.format(p[XKOO],p[XPKOO],p[YKOO],p[YPKOO],p[ZKOO],p[ZPKOO],p[EKOO],p[SKOO])
        return s
    # soll = None

# default track-point                 x   x'  y   y'  z   z'           Tk                 1   s   1
# SollTrack = Track(start=np.array([ 0., 0., 0., 0., 0., 0., PARAMS['sollteilchen'].tkin, 1., 0., 1.]))

def print_progress_bar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '='):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total:
        print()
#
# Sample Usage of print_progress_bar(...)
#
#
# from time import sleep
#
# # make a list
# items = list(range(0, 57))
# i = 0
# l = len(items)
#
# # Initial call to print 0% progress
# print_progress_bar(i, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
# for item in items:
#     # Do stuff...
#     sleep(0.1)
#     # Update Progress Bar
#     i += 1
#     print_progress_bar(i, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
#
# # Sample Output
# Progress: |█████████████████████████████████████████████-----| 90.0% Complete

def track_soll(lattice):
    """
    Tracks the reference particle through the lattice and redefines the lattice element parameters to
    adapted to the energy of the accellerated reference particle.
    """
    # sollteilchen track-point          x   x'  y   y'  z   z'          Tk                   1   s   1
    soll_track = Track(start=np.array([ 0., 0., 0., 0., 0., 0., PARAMS['sollteilchen'].tkin, 1., 0., 1.]))
      #track of reference particle
    for ipos in lattice.seq:
        element,s0,s1 = ipos
        # DEBUG_MODULE('\n{}\t(#{}, pos {:.4f}) label \'{}\''.format(element.__class__,id(element),s0,element.label))
        ti = soll_track.last()                #track: at entrance
        # DEBUG_MODULE('\t\ti >>',Track.string(ti))
        element.adjust_energy(ti[EKOO])    #enery adaptation
        tf = element.soll_map(ti)             #track: at exit
        # DEBUG_MODULE('\t\tf >>',Track.string(tf))
        soll_track.append(tf)
        # deltaE = tf[EKOO] - ti[EKOO]
        # DEBUG_MODULE('\t\tf >>',Track.string(tf),' deltaE[KeV] >>',deltaE*1.e3)
    # DEBUG_MODULE('complete track\n{}'.format(soll_track.points_string()))
    # DEBUG_MODULE('soll track(i)\n{}'.format(soll_track.first_str()))
    # DEBUG_MODULE('soll track(f)\n{}'.format( soll_track.last_str()))
    return soll_track

def track(lattice,bunch):
    """
    Tracks a bunch of particles through the lattice
    lattice: a list of elements a.k.a. _matrix'es
    bunch: a list of independent Tracks
    """
    from time import sleep
    print_progress_bar(0, bunch.nb_tracks(), prefix = 'Progress:', suffix = 'Complete', length = 50)
    for (count,particle_track) in enumerate(bunch.tracks()):
        ti = particle_track.first()
        for ipos in lattice.seq:
            element,s0,s1 = ipos
            tf = element.matrix.dot(ti)      #track through!
            if isinstance(element,ELM.MRK):
                # if count == 1: DEBUG_MODULE('tf in track({}) >>'.format(count)+Track.string(tf))
                particle_track.append(tf)
                # deltaE = tf[EKOO] - ti[EKOO]
                # DEBUG_MODULE('\t\tf >>',Track.string(tf),' deltaE[KeV] >>',deltaE*1.e3)
            ti = tf
        sleep(1.0e-3)
        print_progress_bar(count, bunch.nb_tracks(), prefix = 'Progress:', suffix = 'Complete', length = 50)
        # DEBUG_MODULE('complete track\n{}'.format(particle_track.points_str()))
        # DEBUG_MODULE('FIRST: {}'.format(particle_track.first_str()))
        # DEBUG_MODULE('{} LAST: {}'.format(count,particle_track.last_str()))
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    print("tracks.py: sorry - nothing todo")
