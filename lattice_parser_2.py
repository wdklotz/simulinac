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
ParserResult = namedtuple('ParserResult','FLAGS, PARAMETERS, ELEMENTS, LATTICE, LAT_ELMIDs, ELMIDs')

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

def unnest(dictionary,scalarp,PARTS,result=None):
    if result is None: result = []
    for k,v in dictionary.items():
        DEBUG_OFF('{}:{} is scalar? {}'.format(k,v,scalarp(v)))
        if not scalarp(v):
            for key in v:
                unnest({key:PARTS[key]},scalarp,PARTS,result)  # recursive call !!
        else:
            result.append(v)
    return result

def remove_duplicates(elementIDs):
    """Remove duplicate elementIDs in a sequence"""
    DEBUG_OFF(elementIDs)
    seen = set()
    for elementID in elementIDs:
        if elementID in seen:
            continue
        else:
            seen.add(elementID)
    return list(seen)

# ------------ parser body    ------------ parser body    ------------ parser body    ------------ parser body    ------------ parser body    
# ------------ parser body    ------------ parser body    ------------ parser body    ------------ parser body    ------------ parser body    
# ------------ parser body    ------------ parser body    ------------ parser body    ------------ parser body    ------------ parser body    
def parse(in_data):
    DEBUG_OFF(in_data)

    PARTS = {}

    DEBUG_OFF(HR)
    DEBUG_OFF('FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS    FLAGS    FLAGS')
    DEBUG_OFF(HR)
    flags = in_data['FLAGS']
    DEBUG_OFF(flags)

    DEBUG_OFF(HR)
    DEBUG_OFF('PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS')
    DEBUG_OFF(HR)
    parameters = in_data['PARAMETERS']
    DEBUG_OFF(parameters)

    DEBUG_OFF(HR)
    DEBUG_OFF('ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS')
    DEBUG_OFF(HR)
    elements = in_data['ELEMENTS']
    DEBUG_OFF(elements)

    DEBUG_OFF(HR)
    DEBUG_OFF('NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  ')
    DEBUG_OFF(HR)
    elementIDs = list(elements)
    DEBUG_OFF(elements)

    DEBUG_OFF(HR)
    DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
    DEBUG_OFF(HR)
    for elementID in elementIDs:
        PARTS[elementID] = elementID
    DEBUG_OFF(PARTS)

    DEBUG_OFF(HR)
    DEBUG_OFF('SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS')
    DEBUG_OFF(HR)
    segments = in_data['SEGMENTS']
    DEBUG_OFF(segments)
    segments=expand(segments)
    DEBUG_OFF(segments)

    DEBUG_OFF(HR)
    DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
    DEBUG_OFF(HR)
    for k,v in segments.items():
        PARTS[k]=v
    DEBUG_OFF(PARTS)

    DEBUG_OFF(HR)
    DEBUG_OFF('CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS')
    DEBUG_OFF(HR)
    cells = in_data['CELLS']
    DEBUG_OFF(cells)
    cells = expand(cells)
    DEBUG_OFF(cells)

    DEBUG_OFF(HR)
    DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
    DEBUG_OFF(HR)
    for k,v in cells.items():
        PARTS[k]=v
    DEBUG_OFF(PARTS)

    DEBUG_OFF(HR)
    DEBUG_OFF('SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS')
    DEBUG_OFF(HR)
    sections = in_data['SECTIONS']
    DEBUG_OFF(sections)     # this shows what the YAML parser has done
    sections = expand(sections)# apply_NTIMES(sections)  # expand 'ITEMS'
    DEBUG_OFF(sections)

    DEBUG_OFF(HR)
    DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
    DEBUG_OFF(HR)
    for k,v in sections.items():
        PARTS[k]=v
    DEBUG_OFF(PARTS)

    DEBUG_OFF(HR)
    DEBUG_OFF('LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE')
    DEBUG_OFF(HR)
    lattice = in_data['LATTICE']
    DEBUG_OFF(lattice)     # this shows what the YAML parser has done
    lattice = expand(lattice)# apply_NTIMES(sections)  # expand 'ITEMS'
    DEBUG_OFF(lattice)

    DEBUG_OFF(HR)
    DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
    DEBUG_OFF(HR)
    for k,v in lattice.items():
        PARTS[k]=v
    DEBUG_OFF(PARTS)

    lat_elmIDs = unnest(lattice, is_string_like, PARTS)
    DEBUG_OFF(lat_elmIDs)

    elmIDs = remove_duplicates(lat_elmIDs)
    DEBUG_OFF(elmIDs)

    # ParserResult = namedtuple('ParserResult','FLAGS, PARAMETERS, ELEMENTS, LATTICE, LAT_ELMIDs, ELMIDs')
    ParserResult.FLAGS          = flags
    ParserResult.PARAMETERS     = parameters
    ParserResult.ELEMENTS       = elements
    ParserResult.LATTICE        = lattice
    ParserResult.LAT_ELMIDs     = lat_elmIDs
    ParserResult.ELMIDs         = elmIDs
    return ParserResult

def test0(input_file):
    print(HR+'> test0')
    with open(input_file, 'r') as f:
        in_data = yaml.load(f,Loader=yaml.Loader)
    results = parse(in_data)
    DEBUG_OFF(results.FLAGS)
    DEBUG_OFF(results.PARAMETERS)
    DEBUG_ON(results.ELEMENTS)
    DEBUG_OFF(results.LATTICE)
    DEBUG_ON(results.LAT_ELMIDs)
    DEBUG_ON(results.ELMIDs)

if __name__ == '__main__':
    # test0('yml/lattice_parser_2testIN.yml')
    test0('yml/lattice_parser2_compatIN.yml')
