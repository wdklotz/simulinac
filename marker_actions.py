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
import matplotlib.pyplot as plt

import elements as ELM
from setutil import PARAMS,Ktp

class PoincareAction(ELM.MRK):
    """
    This marker-action will be used by the 'scatter' action
    """
    def __init__(self, label='PSC', prefix='', particle=PARAMS['sollteilchen'], position=[0, 0, 0], aperture=None, next=None, prev=None, abszisse='z', ordinate='zp'):
        super().__init__(label=label, particle=particle, position=position, aperture=aperture, next=next, prev=prev, action='scatter')
        # all points for this scatter-marker
        self.tpoints  = []
        self.prefix   = prefix
        krows         = dict(x=Ktp.x.value, xp=Ktp.xp.value, y=Ktp.y.value, yp=Ktp.yp.value, z=Ktp.z.value, zp=Ktp.zp.value)
        self.abszisse = abszisse
        self.ordinate = ordinate
        self.xaxis    = krows[abszisse]
        self.yaxis    = krows[ordinate]

    def add_track_point(self,track_point):
        self.tpoints.append(track_point)
        
    def do_action(self,number,xmax,ymax):
        # box is text in upper left in axes coords (remark: these are matplotlib.patch.Patch properties)
        symbols = ('x', "x'", 'y', "y'", 'z', "$\Delta$p/p")
        box = '{}-{}'.format(symbols[self.xaxis],symbols[self.yaxis])
        box = '{0}, {2:.2f}[m], {1}p'.format(box,len(self.tpoints),self.position[1])
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        ax = plt.subplot(label=number)
        ax.text(0.05, 0.95, box, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)
        plt.xlim([-xmax,xmax])
        plt.ylim([-ymax,ymax])
        ax.autoscale(enable=False,axis='both')
        x = [tp()[self.xaxis] for tp in self.tpoints]
        y = [tp()[self.yaxis] for tp in self.tpoints]
        ax.scatter(x,y,s=1)
        plt.savefig('{}/poincare_section_{:03d}.png'.format(self.prefix,number))

    def adjust_energy(self, tkin):
        _params = self._params
        self.__init__(label=self.label, prefix=self.prefix, particle=self.particle(tkin), position=self.position, aperture=self.aperture, next=self.next, prev=self.prev, abszisse=self.abszisse, ordinate=self.ordinate)
        self._params = _params
        return self
