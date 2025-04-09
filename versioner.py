#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
__version__='v11.0.3'
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
import sys,os
import glob
import re

""" Utility to change the __version__ global in all *.py-files to a new value """

def replace(fIN):
    fIN.seek(0)
    pattern = "__version__=[\"\'v]\d{1,2}[.]\d{1,2}[.]\S{1,5}[\"\']"
    lines = fIN.readlines()
    for cnt,line in enumerate(lines):
        match = re.search(pattern,line)
        if match:
            print(f'match is {match.group()} in file {fIN.name}')
            lines[cnt] = "__version__='{}'\n".format(new_version)

    with open('new_versions/'+fIN.name,'w',encoding="utf-8") as fOUT:
        fOUT.writelines(lines)
    fOUT.close()

def main():
    py_files = [f for f in glob.glob('*.py')]
    for cnt,py_file in enumerate(py_files):
        # if cnt != 0 : continue
        with open(py_file,"r",encoding="utf-8") as fIN:
            for i in range(10):
                line = fIN.readline()
                if line.find('__version__') != -1:    replace(fIN)  # copy file and replace line
            fIN.close()
    print("DONE: new version {} ==> {} *.py-files stored in directory ./new_versions".format(new_version,cnt))

if __name__ == '__main__':
    new_version = sys.argv[1] if len(sys.argv) == 2 else None
    if new_version == None:
        print("usage: versioner.py v<new version>")
        sys.exit(1)
    main()



