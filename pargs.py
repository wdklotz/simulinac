##!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
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

    You should have received a copy of the GNU General Public Licensedir
    along with SIMULINAC.  If not, see <http://www.gnu.org/licenses/>.
"""
import sys

def input_error():
    print("Input Error!")
    print("Usage: python (simu.py | tracker.py) [(file | --tmpl number --run number)]")
    sys.exit(1)

def pargs(args):
    # print(args)
    proc    = args[0]   
    # either simu.py or tracker.py
    if proc != 'simu.py' and proc != 'tracker.py':
        input_error()
    command = 'python'
    file    = 'yml/simuIN.yml' if proc == 'simu.py' else 'yml/trackIN.yml'
    # args normalized
    Args = {  'mode'  : 'no_m4',
              'proc'  : proc,
              'file'  : file,
              'tmpl'  : '',
              'macro' : ''
            }
    if len(args) == 1:
        # python simu.py
       pass
    elif len(args) == 2:
        # python simu.py file
        file = args[1]
        Args['file'] = file
    elif len(args) == 3:
        # python simu.py --tmpl Number
        opt  = args[1]
        tmpl = args[2]
        Args['mode'] = 'm4'
        if opt == '--tmpl':
            Args['tmpl'] = 'yml/tmpl_{}.yml'.format(tmpl)
            Args['macro']= 'yml/macros_{}.sh'.format(tmpl)
        else: 
            input_error()
    elif len(args) == 5:
        # python simu.py --tmpl Number --run Number
        opt1 = args[1]
        num1 = args[2]
        opt2 = args[3]
        num2 = args[4]
        Args['mode'] = 'm4'
        if   opt1 == '--tmpl' and opt2 == '--run':
            Args['tmpl']   = 'yml/tmpl_{}.yml'.format(num1)
            Args['macro']  = 'yml/macros_{}.{}.sh'.format(num1,num2)
        elif opt1 == '--run' and opt2 == '--tmpl':
            Args['tmpl']   = 'yml/tmpl_{}.yml'.format(num2)
            Args['macro']  = 'yml/macros_{}.{}.sh'.format(num2,num1)
        else:
            input_error()
    else:
        input_error()

    return Args