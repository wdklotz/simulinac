from elements import MRK
from setutil import Twiss, DEBUG_ON, DEBUG_OFF, elli_sxy_action

class PSmarker(MRK):
    def __init__(self, label="Mark1", position=[1,1,1], values=(0.5,0.5),which_action="transvers"):
        super().__init__(label=label, position=position)
        self.label = label
        self.position = position
        self.twiss_values = values
        self.which = which_action
        self.my_actions=dict(transvers=self.action_transvers, no_plot=self.action_no_plot)  # dispatcher

    def do_action(self):
        self.my_actions.get(self.which)()   # dispatch which_action

    def action_transvers(self):
            elli_sxy_action(None,on_injection=True)
            
    def action_no_plot(self):
            DEBUG_ON(self.__dict__)

def test0():
    print('----------------------------------------- test0')
    action0 = PSmarker(which_action='no_plot')
    action1 = PSmarker(label='MyMarker', which_action='no_plot')
    markr = MRK(action=action0)
    markr.add(action1)
    markr.do_actions()
def test1():
    print('----------------------------------------- test1')
    action0 = PSmarker()
    markr = MRK(action=action0)
    markr.do_actions()
    plt.show()
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    test0()
    test1()
