import sys
import yaml
from collections import namedtuple
import pprint, inspect

def PRINT_PRETTY(obj=None):
    file = inspect.stack()[0].filename
    if obj != None: print('DEBUG_ON ==============>  '+file)
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj=None):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

HR = '============================================================================================='

ParserResult = namedtuple('ParserResult','DESCRIPTOR, FLAGS, PARAMETERS, ELEMENTS, LATTICE, LAT_ELMIDs, ELMIDs')

def is_string_like(obj):
    try: obj+'x'
    except TypeError: return False
    else: return True
def nlists(segments):
    """ make n*[...] """
    _segments = {}
    for k,v in segments.items():
        n = int(v[0])
        l = v[1]
        _segments[k] = n*l
    return _segments
def unnest(dictionary,scalarp,PARTS,result=None):
    """ unnest tree of lists down to leaves  (recursive function)"""
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
__parserResult = None    # parser results kept as global module variable
def parse(in_data=None):
    """
    # Module with global.
    # Simpler than the Singleton pattern or BORG idiom.
    # Recommended by Python Cookbook pp.209
    """
    global __parserResult
    if in_data != None:
        DEBUG_OFF(in_data)
        PARTS = {}            
        # descriptor = in_data['DESCRIPTOR'] if 'DESCRIPTOR' in in_data else None
        descriptor = in_data.get('DESCRIPTOR')  # doing the dame as line above
        DEBUG_OFF(descriptor)

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
        segments=nlists(segments)
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
        cells = nlists(cells)
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
        sections = nlists(sections)# apply_NTIMES(sections)  # nlists 'ITEMS'
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
        lattice = nlists(lattice)
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

        """
        DESCRIPTOR: some text to describe the INput
        FLAGS: dict of all flags
        PARAMETERS: dict of all parameters
        ELEMENTS: dict of all ELEMENT types
        LATTICE: list result of appling n*[] for lattice
        LAT_ELMIDs: list result of all elementIDs after unnesting
        ELMIDs: list result of element IDs after removal of duplicate items in LAT_ELMIDs
        __parserResult.DESCRIPTOR     = descriptor
        __parserResult.FLAGS          = flags
        __parserResult.PARAMETERS     = parameters
        __parserResult.ELEMENTS       = elements
        __parserResult.LATTICE        = lattice
        __parserResult.LAT_ELMIDs     = lat_elmIDs   
        __parserResult.ELMIDs         = elmIDs   #TODO not neede anymore?
        """
        __parserResult = ParserResult(descriptor,flags,parameters,elements,lattice,lat_elmIDs,elmIDs)
    else:
        pass
    return __parserResult

def test0(input_file):
    print(HR+'> test0')
    with open(input_file, 'r') as f:
        in_data = yaml.load(f,Loader=yaml.Loader)
    results = parse(in_data)    # this parses the IN-file and returns the  results of parsing the IN-file
    results = parse()           # this only returns the results of parsing the IN-file
    DEBUG_ON(results.DESCRIPTOR)
    DEBUG_OFF(results.FLAGS)
    DEBUG_OFF(results.PARAMETERS)
    DEBUG_OFF(results.ELEMENTS)
    DEBUG_OFF(results.LATTICE)
    DEBUG_OFF(results.LAT_ELMIDs)
    DEBUG_OFF(results.ELMIDs)

if __name__ == '__main__':
    args = sys.argv
    input_file = args[1]
    test0(input_file)
