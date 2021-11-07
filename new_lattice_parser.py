import sys
import yaml
from collections import namedtuple
import pprint

def PRINT_PRETTY(obj):
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj):
    pass
DEBUG_ON  = PRINT_PRETTY
DEBUG_OFF = PASS

HR = '============================================================================================='

def apply_NTIMES(segs):
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

def unnest_ITEMS(lattice,elements):
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
            DEBUG_OFF(itm)
            unnest_ITEMS(itm,elements)
    # is next item a hash?
    elif isinstance(root,dict):
        # print('dict_len: {}'.format(len(root)))
        for key in list(root):
            if key == "DESC": continue
            # print('key: {}'.format(key))
            DEBUG_OFF(root[key])
            unnest_ITEMS(root[key],elements)
    # is next item a string?
    # elif isinstance(root,str): works, but Cookbook prefers
    elif is_string_like(root):
        # print('element {}'.format(root))
        elements.append(root)
    else:
        return

def remove_duplicates(sectionIDs,sections):
    """Remove duplicate elementIDs in a section"""
    # results = namedtuple('InputParseResult',['SECTIONS','LATTICE','FLAGS','PARAMETERS','ELEMENTS'])
    # results = parse(PARAMS['in_data'])
    # sections = results.SECTIONS
    DEBUG_OFF(sections)
    eIDsps = {}
    for sctn in sectionIDs:
        DEBUG_OFF(sctn)
        elementIDs = sections.get(sctn)
        seen = set()
        for element in elementIDs:
            if element in seen:
                continue
            else:
                seen.add(element)
        eIDsps[sctn] = seen
    DEBUG_OFF(eIDsps)
    return eIDsps

def parse(in_data):
    DEBUG_OFF(in_data)
    # print(HR)
    # print('FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS    FLAGS    FLAGS')
    # print(HR)
    flags = in_data['FLAGS']
    DEBUG_OFF(flags)
    # print(HR)
    # print('PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS')
    # print(HR)
    parameters = in_data['PARAMETERS']
    DEBUG_OFF(parameters)
    # print(HR)
    # print('ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS')
    # print(HR)
    elements = in_data['ELEMENTS']
    DEBUG_OFF(elements)
    # print(HR)
    # print('NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  ')
    # print(HR)
    nodes = in_data['NODES']
    DEBUG_OFF(nodes)
    # print(HR)
    # print('SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS')
    # print(HR)
    segments = in_data['SEGMENTS']
    DEBUG_OFF(segments)
    apply_NTIMES(segments)   # expand 'ITEMS'
    DEBUG_OFF(segments)
    # print(HR)
    # print('CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  ')
    # print(HR)
    cells = in_data['CELLS']
    DEBUG_OFF(cells)
    # print(HR)
    # print('SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS')
    # print(HR)
    sections = in_data['SECTIONS']
    DEBUG_OFF(sections)     # this shows what the YAML parser has done
    apply_NTIMES(sections)  # expand 'ITEMS'
    DEBUG_OFF(sections)
    
    dict_sections = {} # {ID:[element-ID,...],...} hash of segments (k,v)=(ID,[element,...])
    for k in list(sections):        # sections is dict
        section = sections.get(k)   # section is deeply nested dict of dict with ITEMS
        element_list_per_section = []
        unnest_ITEMS(section,element_list_per_section)  # recursive unnesting ITEMS
        DEBUG_OFF('Number of elements in {}: {}'.format(section.get('DESC'),len(element_list_per_section)))
        DEBUG_OFF(element_list_per_section)
        dict_sections.setdefault(k,element_list_per_section)

    dict_sections_unique = remove_duplicates(list(sections),dict_sections)
    DEBUG_OFF(dict_sections_unique)
    # print(HR)
    # print('LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE')
    # print(HR)
    lattice = in_data['LATTICE']
    DEBUG_OFF(lattice)
    
    ParseResult = namedtuple('InputParseResult',['SECTIONSx','SECTIONSu','LATTICE','FLAGS','PARAMETERS','ELEMENTS'])
    ParseResult.FLAGS      = flags
    ParseResult.PARAMETERS = parameters
    ParseResult.ELEMENTS   = elements
    ParseResult.SECTIONSx  = dict_sections
    ParseResult.SECTIONSu  = dict_sections_unique
    ParseResult.LATTICE    = lattice
    return ParseResult

def test0(input_file):
    print(HR+'> test0')
    with open(input_file, 'r') as f:
        in_data = yaml.load(f,Loader=yaml.Loader)
    # segments, lattice = parse('yml/new-yaml-template.yml')
    results = parse(in_data)
    DEBUG_ON('==================================dict of FLAGS')
    DEBUG_ON(results.FLAGS)
    DEBUG_ON('==================================dict of PARAMETERS')
    DEBUG_ON(results.PARAMETERS)
    DEBUG_ON('==================================dict of ELEMENTS')
    DEBUG_ON(results.ELEMENTS)
    DEBUG_ON('==================================dict of SECTIONS: ITEMS NTIMES expanded in-order')
    DEBUG_ON(results.SECTIONSx)
    DEBUG_ON('==================================dict of SECTIONS with unique ITEMS only')
    DEBUG_ON(results.SECTIONSu)
    DEBUG_ON('==================================LATTICE: dict of actual SECTIONS')
    DEBUG_ON(results.LATTICE)

if __name__ == '__main__':
    test0('yml/tmpl_25.10.2021_new.yml')
