import matplotlib.pyplot as plt
import numpy as np

def NGauss(x,sig,mu):    # Gauss Normalverteilung
    res = np.exp(-(((x-mu)/sig)**2/2.))
    return res

def pre_plt():
    ax  = plt.subplot(111)
    ax.set_title(input_file)
    ax.set_ylabel('Ez0 [MV/m]')
    ax.set_xlabel('z [cm]')
    return ax

def post_plt(ax):
    plt.legend(loc='lower right',fontsize='x-small')
    plt.show()

def display(table,legend):
    zp   = [+float(x[0]) for x in table]
    Ezp  = [+float(x[2]) for x in table]
    zn   = [-float(x[0]) for x in reversed(table)]
    Ezn  = [+float(x[2]) for x in reversed(table)]
    plt.plot(zn+zp,Ezn+Ezp,label=legend)

def test0():
    '''
    Superfish data
    '''
    Ez0_tab = []
    Ez_max = 1.
    with open(input_file,'r') as f:
        lines = list(f)
        for cnt,line in enumerate(lines[41:-2]):
            # print(line,end='')
            stripped    = line.strip()
            (z,sep,aft) = stripped.partition(' ')
            stripped    = aft.strip()
            (R,sep,aft) = stripped.partition(' ')
            stripped     = aft.strip()
            (Ez,sep,aft) = stripped.partition(' ')
            if cnt == 0:
                Ez_max = float(Ez)      # normalization factor
            Ez = float(Ez)/Ez_max
            Ez0_tab.append((z,R,Ez))
    
    # for x in Ez0_tab:
    #     print('z {}[cm]\tR {}[cm]\tEz(z,R) {}[MV/m]'.format(x[0],x[1],x[2]))
    display(Ez0_tab,'sf')

def test1():
    '''
    Gauss'che Normalverteilung
    '''
    z = np.arange(0.,4.,0.04)
    sigm = 1.14
    Ez0_tab = [(x,0.,NGauss(x,sigm,0.)) for x in z]
    cavlen = 2.5
    print('sigma {}[cm], half-cavity length = {}[cm] = {}sigma'.format(sigm,cavlen,cavlen/sigm))

    # for x in Ez0_tab:
    #     print('z {}[cm]\tR {}[cm]\tEz(z,R) {}[MV/m]'.format(x[0],x[1],x[2]))
    display(Ez0_tab,'NG')

if __name__ == '__main__':
    input_file = 'SF_WDK2g44.TBL'
    ax = pre_plt()
    test0()
    test1()
    post_plt(ax)
        