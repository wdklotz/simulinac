import yaml
import warnings
import sys

#---- recursive version
def flatten(lis):
    """Given a list, possibly nested to any level, return it flattened."""
    new_lis = []
    for item in lis:
        if type(item) == type([]):
            new_lis.extend(flatten(item))
        else:
            new_lis.append(item)
    return new_lis

##--------- MAIN
if __name__ == '__main__':
    
    input_file = 'nwlat.yml'
    # input_file = 'learnyaml.yaml'
    with open(input_file,'r') as fileobject:
        try:
            indat = yaml.load(fileobject)
        except Exception as ex:
            warnings.showwarning(
                    'InputError: {} - STOP'.format(str(ex)),
                    UserWarning,
                    'lattice_generator.py',
                    'factory()',
                    )
            sys.exit(1)
    fileobject.close()
    ky = indat.keys()
    for k in ky:
        print()
        print(k)
        klist = indat[k]
        print(klist)
        nlist = flatten(klist)
        if k == 'LATTICE':
            N = nlist[0]
            plist = nlist[1:]
            qlist = plist.copy()
            for i in range(N-1):
                qlist += plist
            nlist=qlist
        print(nlist)
    pass
