
class RefEdge:

	def __init__(self, n1 = 0, n2 = 0, x = 0.0, z = 0.0, p = 0.0):

		self.variance = 0.02

		self.origin = n1
		self.target = n2

		self.edgeX = x
		self.edgeZ = z
		self.edgeP = p

		n1.addEdge(self)
		n2.addEdge(self)

	def getRefConstraint(self):
		return self.edgeX, self.edgeZ, self.edgeP

	def setVariance(self, var):
		self.variance = var



