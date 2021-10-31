import sys
import yaml
import pprint

HR = '============================================================================================='

with open('yml/new-yaml-template.yml', 'r') as f:
    in_data = yaml.load(f,Loader=yaml.Loader)
pp = pprint.PrettyPrinter(indent=2)

# pp.pprint(doc)

# print(HR)
# print('FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS    FLAGS    FLAGS')
# print(HR)
flags = in_data['FLAGS']
# pp.pprint(flags), print()

# print(HR)
# print('PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS')
# print(HR)
parameters = in_data['PARAMETERS']
# pp.pprint(parameters), print()

# print(HR)
# print('ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS')
# print(HR)
elements = in_data['ELEMENTS']
# pp.pprint(elements), print()

# print(HR)
# print('NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  ')
# print(HR)
nodes = in_data['NODES']
# pp.pprint(nodes), print()

def expandItems(segs):
    '''
    An algorithm to expand the 'ITEMS' lists in 'SEGMENTS' and 'SECTIONS' 'NTIMES' times.
    '''
    for k,v in segs.items():
        try:
            items = v['ITEMS']
        except KeyError:
            print('No ITEMS!')    # no ITEMS
            sys.exit(1)
        try:
            nitems = v['NITEMS']    # expand  to get NTIMES 'ITEMS'
            if nitems == 0:         # ignore 0 ITEMS
                pass
            else:                   # assemble NITEMS ITEMS
                itemsNtimes = []
                for _ in range(nitems):
                    itemsNtimes += items
                items = itemsNtimes
            del(segs[k]['NITEMS'])   # we don't need this key anymore
        except KeyError:
            pass
        segs[k]['ITEMS'] = items     # thios is the expanded structure

def extract_elements(lattice,elements):
    '''
    A recursive algorithm to pull the element-IDs out of the structure that
    the YAML parser has produced.
    IN: lattice - the structire form the YAML parser
    OUT: elements - a linear list of elements in the lattice
    '''
    root = lattice
    # is next item a list?
    if isinstance(root,list):
        # print('list_len: {}'.format(len(root)))
        for count, itm in enumerate(root):
            # print('list_index: {}'.format(count))
            # pp.pprint(itm),print()
            extract_elements(itm,elements)
    # is next item a hash?
    elif isinstance(root,dict):
        # print('dict_len: {}'.format(len(root)))
        for key in list(root):
            if key == "DESC": continue
            # print('key: {}'.format(key))
            # pp.pprint(root[key]),print()
            extract_elements(root[key],elements)
    # is next item a string?
    elif isinstance(root,str):
        # print('element {}'.format(root))
        elements.append(root)
    else:
        return

# print(HR)
# print('SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS')
# print(HR)
segments = in_data['SEGMENTS']
# pp.pprint(segments), print()
expandItems(segments)
# pp.pprint(segments), print()

# print(HR)
# print('CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  ')
# print(HR)
cells = in_data['CELLS']
# pp.pprint(cells), print()

# print(HR)
# print('SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS')
# print(HR)
sections = in_data['SECTIONS']
# pp.pprint(sections), print()
expandItems(sections)
# pp.pprint(sections), print()

print(HR)
print('LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE')
print(HR)
lattice = in_data['LATTICE']
# the linear sequence of elements in the lattice as a python list:
listOfElementsinLattice = []
extract_elements(lattice,listOfElementsinLattice)
print('Number of elements in lattice: {}'.format(len(listOfElementsinLattice)))
pprint.PrettyPrinter(compact=True).pprint(listOfElementsinLattice)

