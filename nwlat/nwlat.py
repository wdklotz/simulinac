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
    ks = indat.keys()
    for k in ks:
        print()
        print(k)
        kl = indat[k]
        print(kl)
        print(flatten(kl))
    pass
