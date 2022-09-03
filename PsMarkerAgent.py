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
from math import degrees, sqrt, atan
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import pprint, inspect

def PRINT_PRETTY(obj):
    file = inspect.stack()[0].filename
    print('DEBUG_ON ==============>  '+file)
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

from elements import MRK
from setutil import Twiss, PARAMS, Ktw

class PsMarkerAgent(object):
    """ 
    Is an agent for the Marker node which performs an action. 
    Action is selectable by the which_action argument. 
    Default action is 'transvers'. 
    """
    def __init__(self, position=[1,1,1], values=(0.5,0.5)):
        self.label        = "ps-mkr-agent"
        self.position     = position
        self.twiss_values = values
        self.parent       = None

    def set_parent(self,obj):
        self.parent = obj
    def do_action(self,*args):
        """ the default action: plot transvers ellipses """
        node = self.parent
        if node == None: 
            ellipse_plot(node,on_injection=True,scale=0.5)
        else:
            ellipse_plot(node,on_injection=False,scale=0.5)

def ellipse_plot(node,on_injection=False,scale=1.):   
    def convert(xy,alfa,beta,emit):
        """ convert twiss parameters to plot parameters """
        gamma = (1.+alfa**2)/beta
        tilt = degrees(0.5*atan(2*alfa/(gamma-beta)))  # see CERN's Formelsammlung
        a = sqrt(emit*beta)
        b = sqrt(emit/beta)
        # return matplot.patches.Ellipse(origin=xy,width=a,height=b,angle=tilt) arguments
        return (xy,a,b,tilt)
    #------ function body ------ function body ------ function body ------ function body ------ function body ------ function body 
    #------ function body ------ function body ------ function body ------ function body ------ function body ------ function body 
    """ display x- and y-phase-space ellipses """
    # TODO for tracker.py:  Plot a confidence ellipse of a two-dimensional dataset ==> https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html#sphx-glr-gallery-statistics-confidence-ellipse-py
    if on_injection:
        s = 0.0
        ax = PARAMS['alfax_i']
        bx = PARAMS['betax_i']
        ay = PARAMS['alfay_i']
        by = PARAMS['betay_i']

    else:
        twiss = node.twiss      # alpha, beta, gamma
        s = node.position[1]    # position

        ax = twiss[Ktw.ax]
        bx = twiss[Ktw.bx]
        ay = twiss[Ktw.ay]
        by = twiss[Ktw.by]

    org = (0,0)
    ellix = convert(org,ax,bx,PARAMS['emitx_i'])  # emittance
    elliy = convert(org,ay,by,PARAMS['emity_i'])

    ellipses = [Ellipse(*ellix,color='blue',fill=False),Ellipse(*elliy,color='red',fill=False)]   # classicla Snyder&Courant ellipse

    fig, ax = plt.subplots()
    fig.suptitle('phase-space {{[m],[rad]}} @ s={:6.2f}[m]'.format(s))
    fig.legend(ellipses,("{x,x'}","{y,y'}"),loc=1)
    fig.suptitle('phase-space {{[m],[rad]}} @ s={:6.2f}[m]'.format(s))
    fig.legend(ellipses,("{x,x'}","{y,y'}"),loc=1)
    
    for e in ellipses:
        # ax.add_artist(e)
        ax.add_patch(e)
        e.set_clip_box(ax.bbox)

    x1 = sqrt(PARAMS['emitx_i']*PARAMS['betax_i'])
    x2 = sqrt(PARAMS['emity_i']*PARAMS['betay_i'])
    xmax = max(x1,x2)
    gammax = (1.+PARAMS['alfax_i']**2)/PARAMS['betax_i']
    gammay = (1.+PARAMS['alfay_i']**2)/PARAMS['betay_i']
    y1 = sqrt(PARAMS['emitx_i']*gammax)
    y2 = sqrt(PARAMS['emity_i']*gammay)
    ymax = max(y1,y2)
    plt.xlim(-xmax*scale, xmax*scale)
    plt.ylim(-ymax*scale, ymax*scale)
    # plt.show()    # do not show now! will be done by simu.py
    return

def test0():
    print('----------------------------------------- test0')
    action0 = PsMarkerAgent(which_action='no_plot')
    action1 = PsMarkerAgent(label='MyMarker', which_action='no_plot')
    markr = MRK(agent=action0)    # new Marker with agent
    markr.add(action1)            # ad a second agent
    markr.do_actions()
def test1():
    print('----------------------------------------- test1')
    action0 = PsMarkerAgent()
    markr = MRK(agent=action0)
    markr.do_actions()
    plt.show()
if __name__ == '__main__':
    test0()
    test1()
