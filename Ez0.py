from matplotlib.pyplot import plot,show,legend,figure,subplot,axis

Ez0_tab = []
input_file = 'SF_WDK2g44.TBL'
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
        Ez0_tab.append((z,R,Ez))

for x in Ez0_tab:
    print('z {}[cm]\tR {}[cm]\tEz(z,R) {}[MV/m]'.format(x[0],x[1],x[2]))

zp   = [+float(x[0]) for x in Ez0_tab]
Ezp  = [+float(x[2]) for x in Ez0_tab]
zn   = [-float(x[0]) for x in reversed(Ez0_tab)]
Ezn  = [+float(x[2]) for x in reversed(Ez0_tab)]

ax  = subplot(111)
ax.set_title(input_file)
ax.set_ylabel('Ez0 [MV/m]')
ax.set_xlabel('z [cm]')
ax.plot(zn+zp,Ezn+Ezp)
show()
        