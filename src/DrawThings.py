import pylab
from math import sin, cos

class DrawThings():
	
	def __init__(self, robotParam):
		self.renderView = 0
		self.robotParam = robotParam
		self.count = 0
		self.isSim = False

	def setSim(self):
		self.isSim = True

	def setRenderView(self, renderView):
		self.renderView = renderView
		
	def plotRobotConfiguration(self, poses, clr = (0.0,0.0,0.0)):
		
		if self.isSim:
			return
		
		segWidth = self.robotParam['segWidth']
		segLength = self.robotParam['segLength']
				
		pylab.clf()
		
		for p in poses:
			xP = []
			yP = []
			
			p1 = [p[0] - 0.5*segWidth*sin(p[2]), p[1] + 0.5*segWidth*cos(p[2])]
			p2 = [p[0] + 0.5*segWidth*sin(p[2]), p[1] - 0.5*segWidth*cos(p[2])]
			p3 = [p[0] + segLength*cos(p[2]) + 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) - 0.5*segWidth*cos(p[2])]
			p4 = [p[0] + segLength*cos(p[2]) - 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) + 0.5*segWidth*cos(p[2])]

			xP = [p4[0],p3[0],p2[0],p1[0],p4[0]]
			yP = [p4[1],p3[1],p2[1],p1[1],p4[1]]
			
			pylab.plot(xP,yP, color = clr)
			
		pylab.xlim(-4,9)
		pylab.ylim(-7,5)
		pylab.savefig("pose_%04u.png" % self.count )
		self.count += 1
	
	def plotPoints(self, pnts, color = (0.0,0.0,0.0)):
		pass
	
	def plotLine(self, pnts, color = (0.0,0.0,0.0)):
		pass
	
	def renderPoints(self, pnts, color = (0.0,0.0,0.0)):
		pass

	def renderLines(self, pnts, color = (0.0,0.0,0.0)):
		pass

	def saveView(self, filename):
		if self.renderView != 0:
			self.renderWindow.writeContentsToFile(filename)
			
