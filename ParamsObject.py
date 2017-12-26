class ParamsObject(object):
	"""
	Class. An object that has key:value parameters
	"""
	def __init__(self):
		self._params = {}
		
	def __getitem__(self,k):
		return self._params[k]
		
	def __setitem__(self,k,v):
		self._params[k] = v

class Aclass(ParamsObject,object):
	pass

class Bclass(Aclass):
	pass
	
## main ----------
if __name__ == '__main__':
	o = ParamsObject()
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
	