import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *
from Spring import *
from Anchor import *
from Transition import Transition
from PathCurveFit import *
from pose import ValueStability



class GuidedPoke(Behavior):

	def __init__(self, probe, mapGraph):
		Behavior.__init__(self, probe)
		
		self.mapGraph = mapGraph

		self.spring = Spring(self.probe)
		self.spring.setTopJoint(13)

		self.obstContact = self.mapGraph.obstCallBack
		self.MAX = 10
		self.MIN = 6

		self.rootNode = 19
		self.topJoint = self.MIN

		self.direction = True

		self.state = -1
		
		self.count = 0

		self.pathCurve = PathCurveFit(self.probe, self.mapGraph)

		self.stabilityX = ValueStability(1e-7, 10)
		self.stabilityY = ValueStability(1e-7, 10)
		
		"""
		Method of Operation
		
		1. curl front snake into an accordion
		2. select a frontier point in the current direction
		3. compute a curve from the origin to the target
		4. feed snake segments to make a global curve fit until we reach the target
		5. attempt to extrapolate the curve further and push the boundary
		6. goto 1
		"""

	def setDirection(self, val):
		#print "setting direction", val
		self.direction = val

		if self.direction:
			pass
			#self.curl.setTopJoint(self.topJoint)
		else:
			pass
			#self.curl.setTopJoint(-self.topJoint + 1 + self.probe.numSegs-2)
		
		#self.curl.setDirection(self.direction)


	def getJointPose(self, jointI):
		
		segLength = self.probe.segLength
		segWidth = self.probe.segWidth
		numJoints = self.probe.numSegs-1
		
		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0
	
		if jointI > self.rootNode:
			#joints = range(self.rootNode, self.probe.numSegs-1)
			joints = range(self.rootNode, jointI)
	
			for i in joints:
				xTotal = xTotal + segLength*cos(totalAngle)
				zTotal = zTotal + segLength*sin(totalAngle)
	
				totalAngle = totalAngle - self.probe.getServo(i+1)
				totalAngle = normalizeAngle(totalAngle)
		
		elif jointI < self.rootNode:
			
			joints = range(jointI,self.rootNode)
			joints.reverse()
	
			for i in joints:
				totalAngle = totalAngle + self.probe.getServo(i+1)
				totalAngle = normalizeAngle(totalAngle)
				
				xTotal = xTotal - segLength*cos(totalAngle)
				zTotal = zTotal - segLength*sin(totalAngle)
		
		elif seg == self.rootNode:
			xTotal = 0.0
			zTotal = 0.0

		else:
			print "rootNode =", self.rootNode, "segI = ", segI
			raise

	
		return [xTotal, zTotal, totalAngle]

	def getPathDist(self, path):
		
		totalDist = 0.0
		
		for i in range(len(path)-1):
			p0 = path[i]
			p1 = path[i+1]
			
			dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
			
			totalDist += dist
		
		return totalDist
		

	def step(self):
		Behavior.step(self)
		
		if not self.isStep():
			return False

		self.resetJoints()

		finalDone = False

		if self.state == -1:
			self.mapGraph.saveLocalMap()
			self.state = 0

		if self.state == 0:

			isDone = self.spring.step()
			joints1 = self.spring.getJoints()
			self.mergeJoints([joints1])

			if isDone:
				self.state = 1

		elif self.state == 1:
			self.localNode = self.mapGraph.getCurrentNode()
			self.frontierMap = self.localNode.getFrontierMap()
			self.navMap = self.localNode.getNavMap()

			self.mapGraph.synch()

			pointNotFound = True
			
			while pointNotFound:
				
				try:
					fPoint = self.frontierMap.selectNextFrontier(self.direction)
				except:	
					self.count = 0
					self.state = 0
					self.mapGraph.saveLocalMap()
					return True
				
				" the root node is our origin "
				if self.direction:
					oPoint = self.getJointPose(13)
				else:
					oPoint = self.getJointPose(25)
					
				
				self.path = self.navMap.computePath(oPoint, fPoint)
				
				maxDist = (13 + 1)*self.probe.segLength;

				pathDist = self.getPathDist(self.path)
				
				print "maxDist =", maxDist
				print "pathDist =", pathDist
				
				if pathDist <= maxDist:
					pointNotFound = False
			
				self.frontierMap.inhibitLocation(fPoint)
				
			self.curve = VoronoiFit(self.path)
			
			self.pathCurve.setCurve(self.curve)
			self.pathCurve.setBoundaries(0,0)
			self.setDirection(True)

			#self.pathCurve.draw()

			self.state = 2
			joints1 = self.spring.getJoints()
			self.mergeJoints([joints1])

			self.stabilityX.reset()
			self.stabilityY.reset()
			self.stabilityX.setThresh(1e-6)
			self.stabilityY.setThresh(1e-6)

		elif self.state == 2:
			self.pathCurve.step()
			self.pathCurve.draw()

			joints1 = self.pathCurve.getJoints()
			joints2 = self.spring.getJoints()
			self.mergeJoints([joints1, joints2])
			
			self.count += 1
			
			if self.direction:
				"boundary, 0 to 13"
				endNode = self.count / 100
				if endNode > 13:
					endNode = 13
				self.pathCurve.setBoundaries(0,endNode)
			else:
				"boundary, 38 to 25"
				startNode = 38 - self.count / 100
				if startNode < 25:
					startNode = 25
				self.pathCurve.setBoundaries(startNode,38)	
			
			if self.count >= 1300:
				if self.direction:
					pose = self.pathCurve.getForeTip()
				else:
					pose = self.pathCurve.getAftTip()
					
				self.stabilityX.addData(pose[0])
				self.stabilityY.addData(pose[1])

				print self.stabilityX.getVar(), self.stabilityY.getVar()

				if self.stabilityX.isStable() and self.stabilityY.isStable():
			
					#if self.count >= 1800:
					self.count = 0
					self.state = 0
					self.obstContact(self.direction)
					self.mapGraph.saveLocalMap()

		joints2 = self.getJoints()

		# joint value
		self.mergeJoints([joints2])


		return finalDone

