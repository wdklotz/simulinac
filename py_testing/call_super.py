class Base():
     def __init__(self,prop):
         self.property = prop 
     def show(self):
        print(self.property)

class Sub(Base):
     def __init__(self,prop):
         super().__init__("base property")
         self.propert = prop
     def show(self):
        print(self.property,end="     ")
        super().show()

sub  = Sub("sub property")
sub.show()


class Parent():
	def show(self):
		print("Inside Parent")
		
class Child(Parent):
	def show(self):
		super().show()
		print("Inside Child")
		
obj = Child()
obj.show()
