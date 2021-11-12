from math import degrees, sqrt, atan
from elements import MRK
from setutil import Twiss, DEBUG_ON, DEBUG_OFF, PARAMS

class PsMarkerAgent(object):
    """ 
    Is an agent for the Marker node which performs an action. 
    Action is selectable by the which_action argument. 
    Default action is 'transvers'. 
    """
    def __init__(self, label="Mark1", position=[1,1,1], values=(0.5,0.5), which_action="transvers"):
        self.label = label
        self.position = position
        self.twiss_values = values
        self.which_action = which_action
        self.my_actions=dict(transvers=self.action_transvers, no_plot=self.action_no_plot)  # dispatcher

    def do_action(self):
        self.my_actions.get(self.which_action)()   # dispatch which_action

    def action_transvers(self):
        """ the default action: plot transvers ellipses """
        ellipse_plot(None,on_injection=True)

    def action_no_plot(self):
        """ action without plotting """
        DEBUG_ON(self.__dict__)

def ellipse_plot(node,on_injection=False):
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
#------ function body ------ function body ------ function body ------ function body ------ function body ------ function body 
    """ display x- and y-phase-space ellipses """
    if on_injection:
        s = 0.0
        ax = PARAMS['alfax_i']
        bx = PARAMS['betax_i']
        ay = PARAMS['alfay_i']
        by = PARAMS['betay_i']

    else:
        twiss,s = node['twiss']

        ax = twiss[Ktw.ax]
        bx = twiss[Ktw.bx]
        ay = twiss[Ktw.ay]
        by = twiss[Ktw.by]

    org = (0,0)
    ellix = convert(org,ax,bx,PARAMS['emitx_i'])
    elliy = convert(org,ay,by,PARAMS['emity_i'])

    ellipses = [Ellipse(*ellix,color='blue',fill=False),Ellipse(*elliy,color='red',fill=False)]

    fig, ax = plt.subplots()
    fig.suptitle('phase-space {{[m],[rad]}} @ s={:6.2f}[m]'.format(s))
    fig.legend(ellipses,("{x,x'}","{y,y'}"),loc=1)

    for e in ellipses:
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)

    x1 = sqrt(PARAMS['emitx_i']*PARAMS['betax_i'])
    x2 = sqrt(PARAMS['emity_i']*PARAMS['betay_i'])
    xmax = max(x1,x2)
    gammax = (1.+PARAMS['alfax_i']**2)/PARAMS['betax_i']
    gammay = (1.+PARAMS['alfay_i']**2)/PARAMS['betay_i']
    y1 = sqrt(PARAMS['emitx_i']*gammax)
    y2 = sqrt(PARAMS['emity_i']*gammay)
    ymax = max(y1,y2)
    # scale = 0.6
    scale = 2.0
    plt.xlim(-xmax*scale, xmax*scale)
    plt.ylim(-ymax*scale, ymax*scale)

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
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    test0()
    test1()
