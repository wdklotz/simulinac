import sys
import yaml
from collections import namedtuple

from setutil import DEB, FLAGS

HR = '============================================================================================='

def expand_single_to_many_items(segs):
    '''
    An algorithm to expand the 'ITEMS' lists in 'SEGMENTS' and 'SECTIONS' 'NTIMES' times.
    '''
    for k,v in segs.items():
        seg = segs[k]   # segment {k:v}
        try:
            items = v['ITEMS']
        except KeyError:
            print('No ITEMS!')      # no ITEMS
            sys.exit(1)
        try:
            nitems = v['NITEMS']    # expand  to get NTIMES 'ITEMS'
            if nitems == 0:         # ignore 0 ITEMS
                pass
            else:                   # repeat ITEMS NITEMS times
                itemsNtimes = []
                for _ in range(nitems):
                    itemsNtimes += items
                items = itemsNtimes
            del(seg['NITEMS'])   # we don't need this key anymore
        except KeyError:
            pass
        seg['ITEMS'] = items     # replace single with expanded structure

def is_string_like(obj):
    try: obj+'x'
    except TypeError: return False
    else: return True

def unnest_elements(lattice,elements):
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
            DEB.get('OFF')(itm)
            unnest_elements(itm,elements)
    # is next item a hash?
    elif isinstance(root,dict):
        # print('dict_len: {}'.format(len(root)))
        for key in list(root):
            if key == "DESC": continue
            # print('key: {}'.format(key))
            DEB.get('OFF')(root[key])
            unnest_elements(root[key],elements)
    # is next item a string?
    # elif isinstance(root,str): works, but Cookbook prefers
    elif is_string_like(root):
        # print('element {}'.format(root))
        elements.append(root)
    else:
        return

def parse(in_data):
    DEB.get('OFF')(in_data)
    # print(HR)
    # print('FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS    FLAGS    FLAGS')
    # print(HR)
    flags = in_data['FLAGS']
    DEB.get('OFF')(flags)
    # print(HR)
    # print('PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS')
    # print(HR)
    parameters = in_data['PARAMETERS']
    DEB.get('OFF')(parameters)
    # print(HR)
    # print('ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS')
    # print(HR)
    elements = in_data['ELEMENTS']
    DEB.get('OFF')(elements)
    # print(HR)
    # print('NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  ')
    # print(HR)
    nodes = in_data['NODES']
    DEB.get('OFF')(nodes)
    # print(HR)
    # print('SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS')
    # print(HR)
    segments = in_data['SEGMENTS']
    DEB.get('OFF')(segments)
    expand_single_to_many_items(segments)   # expand 'ITEMS'
    DEB.get('OFF')(segments)
    # print(HR)
    # print('CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  ')
    # print(HR)
    cells = in_data['CELLS']
    DEB.get('OFF')(cells)
    # print(HR)
    # print('SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS')
    # print(HR)
    sections = in_data['SECTIONS']
    DEB.get('OFF')(sections)
    expand_single_to_many_items(sections)    # expand 'ITEMS'
    DEB.get('OFF')(sections)

    hash_of_sections = {}    # {ID:[element-ID,...],...} hash of segments (k,v)=(ID,[element,...])
    for k in list(sections):        # sections is hash
        section = sections.get(k)   # section is list
        element_list_per_section = []
        unnest_elements(section,element_list_per_section)
        DEB.get('OFF')('Number of elements in {}: {}'.format(section.get('DESC'),len(element_list_per_section)))
        DEB.get('OFF')(element_list_per_section)
        hash_of_sections.setdefault(k,element_list_per_section)

    # print(HR)
    # print('LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE')
    # print(HR)
    lattice = in_data['LATTICE']
    DEB.get('OFF')(lattice)
    
    ParseResult = namedtuple('InputParseResult',['SECTIONS','LATTICE','FLAGS','PARAMETERS','ELEMENTS'])
    ParseResult.FLAGS      = flags
    ParseResult.PARAMETERS = parameters
    ParseResult.ELEMENTS   = elements
    ParseResult.SECTIONS   = hash_of_sections
    ParseResult.LATTICE    = lattice
    return ParseResult

def test0(input_file):
    print(HR+'> test0')
    with open(input_file, 'r') as f:
        in_data = yaml.load(f,Loader=yaml.Loader)
    # segments, lattice = parse('yml/new-yaml-template.yml')
    results = parse(in_data)
    print('\tFLAGS')
    DEB.get('ON')(results.FLAGS)
    print('\tPARAMETERS')
    DEB.get('ON')(results.PARAMETERS)
    print('\tELEMENTS')
    DEB.get('ON')(results.ELEMENTS)
    print('\tSECTIONS')
    DEB.get('ON')(results.SECTIONS)
    print('\tLATTICE')
    DEB.get('ON')(results.LATTICE)

if __name__ == '__main__':
    test0('yml/tmpl_25.10.2021_new.yml')
