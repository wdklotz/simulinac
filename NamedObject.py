class NamedObject:
	"""
	Class. An object that has a name and a label
	"""
	def __init__(self):
		"""
		Object that has tags: name, label, section ...
		"""
		self._label   = "no label"
		self._name    = "no name"
		self._section = "no section"
	
	@property
	def label(self):
		return self._label
		
	@label.setter
	def label(self, value):
		self._label = value

	@property
	def name(self):
		return self._name
		
	@name.setter
	def name(self, value):
		self._name = value
	
	@property
	def section(self):
		return self._section
		
	@name.setter
	def section(self, value):
		self._section = value
	
	
