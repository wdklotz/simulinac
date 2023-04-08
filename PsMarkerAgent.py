#!/Users/klotz/SIMULINAC_env/bin/python
# -*- coding: utf-8 -*-
__version__='v10.23.3'
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
from math import sqrt
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import numpy as np

from elements import MRK
from setutil import Twiss, PARAMS, Ktw, FLAGS, Proton, DEBUG_ON, DEBUG_OFF

class PsMarkerAgent(MRK):
    """ 
    Is an agent for the Marker node which performs 
    a phase-space ellipse plot at the marker's position.
    """
    def __init__(self, label, active, viseo, particle=Proton(PARAMS['injection_energy']), position=[0,0,0], twiss_values=(0.5,0.5)):
        super().__init__(label, active, viseo, particle, position,)
        self.twiss_values = twiss_values
        self.do_action    = self.action if self.active else self.no_action # toggle
    def action(self,*args):
        """ the action: plot transvers ellipses """
        ellipse_plot(self,scale=0.5)

def ellipse_plot(node,scale=1.):   
    def ellipse_and_tranformation(a,b,g,eps,color):
        # https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html
        """
        IN
            a twiss alfa
            b twiss beta
            g twis gamma
            eps emittance
            color
        
        OUT
            sigx sigama abcissa
            sigy sigma ordinate
            ellipse
            transformation
        """
        # covariance matrix from twiss parameters
        cov = np.array([[b,-a],[-a,g]])*eps
        sigx = sqrt(cov[0,0])
        sigy = sqrt(cov[1,1])
        # Pearson correlation coefficient
        pearson = cov[0, 1]/(sigx*sigy)
        ra  = sqrt(1 + pearson)
        rb  = sqrt(1 - pearson)
        # ellipse
        elli = Ellipse((0,0),width=2*rb, height=2*ra, color=color, fill=False)
        # transformation: scaling 1 sigma, rotation 45 deg
        transformation = transforms.Affine2D().rotate_deg(45.).scale(sigx,sigy).translate(0.,0.)
        return (sigx,sigy,elli,transformation)

    #------ function body ------ function body ------ function body ------ function body ------ function body ------ function body 
    """ display x- and y-phase-space ellipses """
    twiss = node.twiss      # alpha, beta, gamma
    s = node.position[1]    # position

    # x,x'
    a = twiss[Ktw.ax]
    b = twiss[Ktw.bx]
    g = twiss[Ktw.gx]
    eps = PARAMS['emitx_i']*1e6
    xlim,xplim,ellix,transx = ellipse_and_tranformation(a,b,g,eps,'blue')
    
    # y,y'
    a = twiss[Ktw.ay]
    b = twiss[Ktw.by]
    g = twiss[Ktw.gy]
    eps = PARAMS['emity_i']*1e6
    ylim,yplim,elliy,transy = ellipse_and_tranformation(a,b,g,eps,'red')

    trafo     = [ transx, transy ]
    ellipses  = [ ellix, elliy ]
    
    fig, ax = plt.subplots(figsize=(5,5))
    fig.suptitle('phase-space {{[mm],[mrad]}} @ s={:6.2f}[m]'.format(s))
    ax.legend(ellipses,("{x,x'}","{y,y'}"),loc='best')
    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)

    for i in range(len(ellipses)):
        e = ellipses[i]
        trans = trafo[i]
        e.set_transform(trans+ax.transData)
        ax.add_patch(e)

    margin = 1.10          # %
    limx=max(xlim,ylim)*margin
    limy=max(xplim,yplim)*margin
    plt.xlim(-limx,limx)
    plt.ylim(-limy,limy)
    # plt.show()    # do not show now! will be done by simu.py
    return

if __name__ == '__main__':
    print('what?')
