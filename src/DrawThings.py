
import ogre.renderer.OGRE as ogre
import pylab
from math import sin, cos

globalID = 0

class DrawThings():
	
	def __init__(self, robotParam):
		self.renderView = 0
		self.robotParam = robotParam
		self.count = 0
		self.isSim = False

		self._allRefNodes = []
		self._allRefEnt = []

		self._allPathNodes = []
		self._allPathEnt = []

	def setSim(self, sceneManager):
		self.sceneManager = sceneManager
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
			
			pylab.plot(xP,yP, linewidth=3, color = clr)
			
		pylab.xlim(-4,9)
		pylab.ylim(-7,5)
		pylab.savefig("pose_%04u.png" % self.count )
		self.count += 1
	
	def plotPoints(self, pnts, color = (0.0,0.0,0.0)):
		pass
	
	def plotLine(self, pnts, color = (0.0,0.0,0.0)):
		pass
	
	def renderPoints(self, pnts, color = (0.0,0.0,0.0)):
		
		global globalID
				
		" delete any existing points that may exist "
		# deference the child nodes now

		for child in self._allRefNodes:
			self.sceneManager.destroySceneNode(child)
	
		for child in self._allRefEnt:
			self.sceneManager.destroyEntity(child)
	
		self._allRefNodes = []
		self._allRefEnt = []
		
		for i in range(len(pnts)):
			## Create the visual reprsentation of active reference nodes
			name = "act_node" + str(globalID)
			entity = self.sceneManager.createEntity(name, "Cube.mesh")
			node = self.sceneManager.getRootSceneNode().createChildSceneNode(name)
			node.attachObject(entity)

			size = ogre.Vector3(0.05,0.05,0.05)
			node.setScale(size)

			pos = pnts[i]
			position = ogre.Vector3(pos[0],0.1,pos[1])
			node.setPosition(position)
			# negative the angle because the bodies are somehow flipped upside down
			node.setOrientation(ogre.Quaternion(-pos[2],ogre.Vector3().UNIT_Y))

			entity.setCastShadows(False)

			entity.setMaterialName("Green")

			#entity.setVisible(False)
			self._allRefNodes.append(node)
			self._allRefEnt.append(entity)
			
			globalID += 1
			
	def drawPath(self, pnts, color = (0.0,0.0,0.0)):
		
		global globalID
				
		" delete any existing points that may exist "
		# deference the child nodes now
		for child in self._allPathNodes:
			self.sceneManager.destroySceneNode(child)
	
		for child in self._allPathEnt:
			self.sceneManager.destroyEntity(child)
	
		self._allPathNodes = []
		self._allPathEnt = []
		
		for i in range(len(pnts)):
			## Create the visual reprsentation of active reference nodes
			name = "act_node" + str(globalID)
			entity = self.sceneManager.createEntity(name, "Cube.mesh")
			node = self.sceneManager.getRootSceneNode().createChildSceneNode(name)
			node.attachObject(entity)

			size = ogre.Vector3(0.05,0.05,0.05)
			node.setScale(size)

			pos = pnts[i]
			position = ogre.Vector3(pos[0],0.1,pos[1])
			node.setPosition(position)

			entity.setCastShadows(False)

			entity.setMaterialName("Red")

			#entity.setVisible(False)
			self._allPathNodes.append(node)
			self._allPathEnt.append(entity)
			
			globalID += 1
				

	def renderLines(self, pnts, color = (0.0,0.0,0.0)):
		pass

	def saveView(self, filename):
		if self.renderView != 0:
			#self.renderView.update()
			self.renderView.writeContentsToFile(filename)
			
