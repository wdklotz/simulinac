#!/Users/klotz/anaconda3/bin/python3.6
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
	
	You should have received a copy of the GNU General Public License
	along with SIMULINAC.  If not, see <http://www.gnu.org/licenses/>.
"""
class DictObject(object):
	"""
	Class. An object that has key:value parameters
	"""
	def __init__(self):
		self._params = {}
		
	def __getitem__(self,k):
		return self._params[k]
		
	def __setitem__(self,k,v):
		self._params[k] = v

class Aclass(DictObject,object):
	pass

class Bclass(Aclass):
	pass
	
## main ----------
if __name__ == '__main__':
	o = DictObject()
	o["name"]="my name is wdk"
	print(o["name"])	

	a = Aclass()
	a["name"]="my name is wdk"
	print(a["name"])
	print(''.join('{}\n'.format(el) for el in Aclass.__mro__))

	b = Bclass()
	b['name']='bbbbbbb'
	print(b['name'])
	print(a['name'] != b['name'])
	print(''.join('{}\n'.format(el) for el in Bclass.__mro__))
	