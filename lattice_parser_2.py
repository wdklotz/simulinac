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

def is_string_like(obj):
    try: obj+'x'
    except TypeError: return False
    else: return True


def expand(segments):
    _segments = {}
    for k,v in segments.items():
        n = int(v[0])
        l = v[1]
        _segments[k] = n*l
    return _segments

def unnest(dictionary,scalarp,result=None):
    if result is None: result = []
    for k,v in dictionary.items():
        DEBUG_OFF('{}:{} is scalar? {}'.format(k,v,scalarp(v)))
        if not scalarp(v):
            for key in v:
                unnest({key:PARTS[key]},scalarp,result)
        else:
            result.append(v)
    return result

PARTS = {}
ELEMENTS = {}
LATTICE = {}

def parse(in_data):
    DEBUG_OFF(in_data)
    global PARTS, ELEMENTS, LATTICE

    DEBUG_ON(HR)
    DEBUG_ON('FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS    FLAGS    FLAGS')
    DEBUG_ON(HR)
    flags = in_data['FLAGS']
    DEBUG_ON(flags)

    DEBUG_ON(HR)
    DEBUG_ON('PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS')
    DEBUG_ON(HR)
    parameters = in_data['PARAMETERS']
    DEBUG_ON(parameters)

    DEBUG_ON(HR)
    DEBUG_ON('ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS')
    DEBUG_ON(HR)
    ELEMENTS = in_data['ELEMENTS']
    DEBUG_ON(ELEMENTS)

    DEBUG_ON(HR)
    DEBUG_ON('NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  ')
    DEBUG_ON(HR)
    nodes = list(ELEMENTS)
    DEBUG_ON(nodes)

    DEBUG_OFF(HR)
    DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
    DEBUG_OFF(HR)
    for node in nodes:
        PARTS[node] = node
    DEBUG_OFF(PARTS)

    DEBUG_OFF(HR)
    DEBUG_OFF('SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS')
    DEBUG_OFF(HR)
    SEGMENTS = in_data['SEGMENTS']
    DEBUG_OFF(SEGMENTS)
    SEGMENTS=expand(SEGMENTS)
    DEBUG_OFF(SEGMENTS)

    DEBUG_OFF(HR)
    DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
    DEBUG_OFF(HR)
    for k,v in SEGMENTS.items():
        PARTS[k]=v
    DEBUG_OFF(PARTS)

    DEBUG_OFF(HR)
    DEBUG_OFF('CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS')
    DEBUG_OFF(HR)
    CELLS = in_data['CELLS']
    DEBUG_OFF(CELLS)
    CELLS = expand(CELLS)
    DEBUG_OFF(CELLS)

    DEBUG_OFF(HR)
    DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
    DEBUG_OFF(HR)
    for k,v in CELLS.items():
        PARTS[k]=v
    DEBUG_OFF(PARTS)

    DEBUG_OFF(HR)
    DEBUG_OFF('SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS')
    DEBUG_OFF(HR)
    SECTIONS = in_data['SECTIONS']
    DEBUG_OFF(SECTIONS)     # this shows what the YAML parser has done
    SECTIONS = expand(SECTIONS)# apply_NTIMES(sections)  # expand 'ITEMS'
    DEBUG_OFF(SECTIONS)

    DEBUG_OFF(HR)
    DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
    DEBUG_OFF(HR)
    for k,v in SECTIONS.items():
        PARTS[k]=v
    DEBUG_OFF(PARTS)

    DEBUG_OFF(HR)
    DEBUG_OFF('LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE')
    DEBUG_OFF(HR)
    LATTICE = in_data['LATTICE']
    DEBUG_OFF(LATTICE)     # this shows what the YAML parser has done
    LATTICE = expand(LATTICE)# apply_NTIMES(sections)  # expand 'ITEMS'
    DEBUG_OFF(LATTICE)

    DEBUG_ON(HR)
    DEBUG_ON('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
    DEBUG_ON(HR)
    for k,v in LATTICE.items():
        PARTS[k]=v
    DEBUG_ON(PARTS)

def test0(input_file):
    print(HR+'> test0')
    with open(input_file, 'r') as f:
        in_data = yaml.load(f,Loader=yaml.Loader)
    results = parse(in_data)
    DEBUG_ON(LATTICE)
    result = unnest(LATTICE,is_string_like)
    DEBUG_ON(result)

if __name__ == '__main__':
    # test0('yml/lattice_parser_2testIN.yml')
    test0('yml/lattice_parser2_compatIN.yml')
