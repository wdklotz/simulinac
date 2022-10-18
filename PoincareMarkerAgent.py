#!/Users/klotz/SIMULINAC_env/bin/python
# -*- coding: utf-8 -*-
__version__='v10.22.5'
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

from elements import MRK
from setutil import PARAMS,Ktp,Proton,FLAGS,DEBUG_ON,DEBUG_OFF

class PoincareMarkerAgent(MRK):
    """ A Marker agent  """
    def __init__(self,label,active,prefix,abscissa,ordinate,particle=Proton(PARAMS['injection_energy']),position=(0.,0.,0.)):
        super().__init__(label, active, particle, position,)
        krows = dict(x=Ktp.x.value, xp=Ktp.xp.value, y=Ktp.y.value, yp=Ktp.yp.value, z=Ktp.z.value, zp=Ktp.zp.value)
        self.tpoints  = []
        self.prefix   = prefix
        self.abscissa = abscissa
        self.ordinate = ordinate
        self.xaxis    = krows[abscissa]
        self.yaxis    = krows[ordinate]
        self.do_action = self.action if FLAGS['maction'] else self.noaction # toggle

    def add_track_point(self,track_point):
        self.tpoints.append(track_point)
    def action(self,*args): 
        """ Generate scatter plot of the poincrecut at this marker and dump it to a file """
        number   = args[0]
        xmax     = args[1]
        ymax     = args[2]
        position = args[3]
        DEBUG_OFF(f'(number,xmax,ymax,position) ({number},{xmax},{ymax},{position})')
        # box is text in upper left in axes coords (remark: these are matplotlib.patch.Patch properties)
        symbols = ('x', "x'", 'y', "y'", 'z', "$\Delta$p/p")
        box = '{}-{}'.format(symbols[self.xaxis],symbols[self.yaxis])
        box = '{0}, {2:.2f}[m], {1}p'.format(box,len(self.tpoints),position[1])
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        ax = plt.subplot(label=number)
        ax.text(0.05, 0.95, box, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)
        plt.xlim([-xmax,xmax])
        plt.ylim([-ymax,ymax])
        ax.autoscale(enable=False,axis='both')
        x = [tp()[self.xaxis] for tp in self.tpoints]
        y = [tp()[self.yaxis] for tp in self.tpoints]
        ax.scatter(x,y,s=1)
        plt.savefig('{}/poincare_cut_{:03d}.png'.format(self.prefix,number))

if __name__ == '__main__':
    print('what?')
    