#-*- coding: utf-8 -*-
import sys,os
import glob
import re

""" Utility to change the __version__ global in all *.py-files to a new value """

def replace(fIN):
    lines = fIN.readlines()
    with open('new_versions/'+fIN.name,'w',encoding="utf-8") as fOUT:
        for line in lines:
            # if line.find('__version__=') != -1:   # the line with the __version__= declaration
            #     line = "__version__='v{}'\n".format(new_version)
            match = rep.match(line)
            if match: 
                # print(match)
                line = "__version__='v{}'\n".format(new_version)
            else:
                pass
            fOUT.write(line)
    fOUT.close()

def main():
    py_files = [f for f in glob.glob('./*.py')]
    f_cnt = 0
    # loop all *.py-files
    for py_file in py_files:
        if py_file.find(sys.argv[0]) == -1:   # skip this script
            with open(py_file,"r",encoding="utf-8") as fIN:
                for i in range(10):
                    line = fIN.readline()
                    if line.find('__version__') != -1:
                        f_cnt += 1
                        fIN.seek(0)
                        replace(fIN)
                fIN.close()
    print("DONE: new version {} ==> {} *.py-files stored in directory ./new_versions".format(new_version,f_cnt))

if __name__ == '__main__':
    # match __version__='vxx.yy.zzz'  xx=main, yy=year, zzz=minor
    rep = re.compile("__version__=['\"]v\d{1,2}[.]\d{1,2}[.]\d{1,3}['\"]") 
    new_version = sys.argv[1] if len(sys.argv) == 2 else None
    if new_version == None:
        print("usage: versioner.py <new version>")
        sys.exit(1)
    main()



