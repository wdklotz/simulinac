import os

from xml_utils.NamedObject   import NamedObject
from xml_utils.TypedObject   import TypedObject

class ParamsDictObject:
	"""
	Class. Object that has a parameters dictionary.
	"""

	def __init__(self):
		"""
		Constructor. Object that has a parameters dictionary.
		"""
		self.__paramsDict = {}

	def addParam(self, key, value):
		"""
		Method. Adds a parameter to the object.
		"""
		self.__paramsDict[key] = value

	def removeParam(self, key):
		"""
		Method. Removes a parameter.
		"""
		if(self.__paramsDict.has_key(key)):
			del  self.__paramsDict[key]

	def setParam(self, key, value):
		"""
		Method. Sets a parameter to the object.
		"""
		self.__paramsDict[key] = value

	def setParamsDict(self, params):
		"""
		Method. Sets an external dictionary as a parameter
		dictionary for this element.
		"""
		self.__paramsDict = params

	def updateParamsDict(self, params):
		"""
		Method. Updates the dictionary with external dictionary data.
		"""
		self.__paramsDict.update(params)

	def getParam(self, key):
		"""
		Method. Returns requested parameters of the object.
		"""
		if(not self.hasParam(key)):
			msg = "The object does not have a parameter for the key you requested!"
			msg = msg + os.linesep
			msg = msg + "method getParam(self, key)"
			msg = msg + os.linesep
			if(isinstance(self,NamedObject) == True):
				msg = msg + "Name of element = " + self.getName()
				msg = msg + os.linesep
			if(isinstance(self,TypedObject) == True):
				msg = msg + "Type of element = " + self.getType()
				msg = msg + os.linesep
			msg = msg + "key = " + str(key)
			print(msg)
		return self.__paramsDict[key]

	def getParamsDict(self):
		"""
		Method. Returns the whole parameters dictionary.
		"""
		return self.__paramsDict

	def hasParam(self, key):
		"""
		Method. Returns True if the object has a parameter
		for this key. Returns False otherwise.
		"""
# 		return self.__paramsDict.has_key(key)
		return key in self.__paramsDict

	def keys(self):
		""" return the list of the keys for the parameters """
		return self.__paramsDict.keys()
