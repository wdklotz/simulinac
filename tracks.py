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

from setutil import PARAMS,DEBUG,tblprnt
from setutil import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
import elements as ELM

## DEBUG MODULE
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_TRACK      = DEBUG_OFF
DEBUG_SOLL_TRACK = DEBUG_OFF

## Track class
class Track(object):
    """ Track is an ordered list of track-points. A track-point is an array of MDIM coordinates."""

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
        return self.track_points.reshape(self.nb_points_per_track,ELM.MDIM)

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
            str += Track.string(p)+'\n'
        return str
    def string(p):   #single point to string
        s = 'x={:.3e} x\'={:.3e} y={:.3e} y\'={:.3e} z={:.3e} z\'={:.3e}  tk={:.5f} s={:.3f} '.format(p[XKOO],p[XPKOO],p[YKOO],p[YPKOO],p[ZKOO],p[ZPKOO],p[EKOO],p[SKOO])
        return s
    def asTable(self):
        tblheadr = ['    x',"    x'",'    y',"    y'",'    z',"    z'",'  tkin','    s']
        tblrows =[]
        for point in self.points():
            tblrow = [
                '{:8.3f}'.format(point[XKOO]),
                '{:8.3f}'.format(point[XPKOO]),
                '{:8.3f}'.format(point[YKOO]),
                '{:8.3f}'.format(point[YPKOO]),
                '{:8.3f}'.format(point[ZKOO]),
                '{:8.3f}'.format(point[ZPKOO]),
                '{:8.3f}'.format(point[EKOO]),
                '{:8.3f}'.format(point[SKOO]),
                ]
            tblrows.append(tblrow)
        return tblprnt(tblheadr,tblrows)

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

def track_soll(lattice):
    """
    Tracks the reference particle through the lattice and redefines the lattice element parameters to
    adapted to the energy of the accellerated reference particle.
    """
    soll_track = Track(start=np.array([ 0., 0., 0., 0., 0., 0., PARAMS['sollteilchen'].tkin, 1., 0., 1.]))
    for element in lattice.seq:
        DEBUG_SOLL_TRACK(element,' pos {:.4f} label "{}"'.format(element.position[1],element.label))
        ti = soll_track.last() #track: at entrance
        DEBUG_SOLL_TRACK('track_soll(i) ',Track.string(ti))
        # DEBUG_SOLL_TRACK('track_soll: complete track\n{}'.format(soll_track.points_str()))
        """ energy adjustment """
        element.adjust_energy(ti[EKOO])
        # DEBUG_SOLL_TRACK('track_soll: complete track\n{}'.format(soll_track.points_str()))
        """ element mapping """
        tf = element.soll_map(ti) #track: at exit
        DEBUG_SOLL_TRACK('track_soll(f) ',Track.string(tf))
        soll_track.append(tf)
        # DEBUG_SOLL_TRACK('track_soll: complete track\n{}'.format(soll_track.points_str()))
    DEBUG_SOLL_TRACK('track_soll: complete track\n{}'.format(soll_track.points_str()))
    DEBUG_SOLL_TRACK('track_soll(first) {}'.format(soll_track.first_str()))
    DEBUG_SOLL_TRACK('track_soll(last)  {}'.format( soll_track.last_str()))
    return soll_track

#todo: @@@track must be rewritten to include non-linear tracking@@@
def track(lattice,bunch):
    """
    Tracks a bunch of particles through the lattice using matrices
    lattice: a list of elements
    bunch: a list of particle tracks
    """
    raise UserWarning("track(lattice,bunch) not up-to-date!!!")
    from time import sleep
    print_progress_bar(0, bunch.nb_tracks(), prefix = 'Progress:', suffix = 'Complete', length = 50)
    for (count,particle_track) in enumerate(bunch.tracks()):
        ti = particle_track.first()
        for element in lattice.seq:
            tf = element.matrix.dot(ti)      # track with matrix!
            if isinstance(element,ELM.MRK):
                if count == 1: 
                    DEBUG_TRACK('track(f)#{}) {}'.format(count,Track.string(tf)))
                particle_track.append(tf)
                DEBUG_TRACK('track(f)#{} {} dE[KeV] {}'.format(count,Track.string(tf),(tf[EKOO] - ti[EKOO])*1.e3))
            ti = tf
        DEBUG_TRACK('track: complete track\n{}'.format(soll_track.points_str()))
        DEBUG_TRACK('track(first) {}'.format(soll_track.first_str()))
        DEBUG_TRACK('track(last)  {}'.format( soll_track.last_str()))

        sleep(1.0e-3)
        print_progress_bar(count, bunch.nb_tracks(), prefix = 'Progress:', suffix = 'Complete', length = 50)
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    print("tracks.py: sorry - nothing todo")
#
# Sample Usage of print_progress_bar(...)
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
