class NamedObject:
	"""
	A class to give objects tags like name a.k.a. label, section ...
	"""
	@property
	def label(self):
		return self._label
		
	@label.setter
	def label(self, value):
		self._label = value

	@property
	def name(self):
		return self._label
		
	@name.setter
	def name(self, value):
		self._label = value
	
	@property
	def section(self):
		return self._section
		
	@section.setter
	def section(self, value):
		self._section = value
	
	
