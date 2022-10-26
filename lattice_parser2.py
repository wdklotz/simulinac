#!/Users/klotz/SIMULINAC_env/bin/python
# -*- coding: utf-8 -*-
__version__='vv10.22.7'
"""
Copyright 2015 Wolf-Dieter Klotz <wdklotz@gmail.com>
This file is part of the SIMULINAC code

    SIMULINAC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    SIMULINAC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SIMULINAC.  If not, see <http://www.gnu.org/licenses/>.
"""
import sys
import yaml
from collections import namedtuple
import setutil as UTIL

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
        UTIL.DEBUG_OFF('{}:{} is scalar? {}'.format(k,v,scalarp(v)))
        if not scalarp(v):
            for key in v:
                unnest({key:PARTS[key]},scalarp,PARTS,result)  # recursive call !!
        else:
            result.append(v)
    return result
def remove_duplicates(elementIDs):
    """Remove duplicate elementIDs in a sequence"""
    UTIL.DEBUG_OFF(elementIDs)
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
        UTIL.DEBUG_OFF(in_data)
        PARTS = {}            
        # descriptor = in_data['DESCRIPTOR'] if 'DESCRIPTOR' in in_data else None
        descriptor = in_data.get('DESCRIPTOR')  # doing the dame as line above
        UTIL.DEBUG_OFF(descriptor)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS  FLAGS    FLAGS    FLAGS')
        UTIL.DEBUG_OFF(HR)
        flags = in_data['FLAGS']
        UTIL.DEBUG_OFF(flags)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS  PARAMETERS')
        UTIL.DEBUG_OFF(HR)
        parameters = in_data['PARAMETERS']
        UTIL.DEBUG_OFF(parameters)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS  ELEMENTS')
        UTIL.DEBUG_OFF(HR)
        elements = in_data['ELEMENTS']
        UTIL.DEBUG_OFF(elements)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  NODES  ')
        UTIL.DEBUG_OFF(HR)
        elementIDs = list(elements)
        UTIL.DEBUG_OFF(elements)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
        UTIL.DEBUG_OFF(HR)
        for elementID in elementIDs:
            PARTS[elementID] = elementID
        UTIL.DEBUG_OFF(PARTS)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS  SEGMENTS')
        UTIL.DEBUG_OFF(HR)
        segments = in_data['SEGMENTS']
        UTIL.DEBUG_OFF(segments)
        segments=nlists(segments)
        UTIL.DEBUG_OFF(segments)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
        UTIL.DEBUG_OFF(HR)
        for k,v in segments.items():
            PARTS[k]=v
        UTIL.DEBUG_OFF(PARTS)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS  CELLS')
        UTIL.DEBUG_OFF(HR)
        cells = in_data['CELLS']
        UTIL.DEBUG_OFF(cells)
        cells = nlists(cells)
        UTIL.DEBUG_OFF(cells)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
        UTIL.DEBUG_OFF(HR)
        for k,v in cells.items():
            PARTS[k]=v
        UTIL.DEBUG_OFF(PARTS)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS  SECTIONS')
        UTIL.DEBUG_OFF(HR)
        sections = in_data['SECTIONS']
        UTIL.DEBUG_OFF(sections)     # this shows what the YAML parser has done
        sections = nlists(sections)# apply_NTIMES(sections)  # nlists 'ITEMS'
        UTIL.DEBUG_OFF(sections)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
        UTIL.DEBUG_OFF(HR)
        for k,v in sections.items():
            PARTS[k]=v
        UTIL.DEBUG_OFF(PARTS)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE  LATTICE')
        UTIL.DEBUG_OFF(HR)
        lattice = in_data['LATTICE']
        UTIL.DEBUG_OFF(lattice)     # this shows what the YAML parser has done
        lattice = nlists(lattice)
        UTIL.DEBUG_OFF(lattice)

        UTIL.DEBUG_OFF(HR)
        UTIL.DEBUG_OFF('PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS  PARTS')
        UTIL.DEBUG_OFF(HR)
        for k,v in lattice.items():
            PARTS[k]=v
        UTIL.DEBUG_OFF(PARTS)

        lat_elmIDs = unnest(lattice, is_string_like, PARTS)
        UTIL.DEBUG_OFF(lat_elmIDs)

        elmIDs = remove_duplicates(lat_elmIDs)
        UTIL.DEBUG_OFF(elmIDs)

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
    UTIL.DEBUG_OFF(results.DESCRIPTOR)
    UTIL.DEBUG_OFF(results.FLAGS)
    UTIL.DEBUG_OFF(results.PARAMETERS)
    UTIL.DEBUG_OFF(results.ELEMENTS)
    UTIL.DEBUG_OFF(results.LATTICE)
    UTIL.DEBUG_OFF(results.LAT_ELMIDs)
    UTIL.DEBUG_ON(results.ELMIDs)
if __name__ == '__main__':
    args = sys.argv
    input_file = args[1] if len(args) >1 else 'unittests/simuIN.yml'
    test0(input_file)
