import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from copy import *
from common import *
from Behavior import *
from FrontAnchorFit import FrontAnchorFit
from FastLocalCurveFit import FastLocalCurveFit
from HoldTransition import HoldTransition
from HoldSlideTransition import HoldSlideTransition
from AdaptiveAnchorCurve import AdaptiveAnchorCurve
from BackConcertinaCurve import BackConcertinaCurve
from GlobalCurveFit import GlobalCurveFit
from GlobalCurveSlide import GlobalCurveSlide

"""
TODO

1.  Check condition in FrontAnchorFit when curve goes beyond joint 0
2.  If joint 0 passed and anchor is loose, start over with higher splice joint
3.  Changing splice joint will shorten tail of BackConcertinaCurve
4.  Add the front anchoring state and then the concertina gait state to behavior control


"""

class PathStep(Behavior):

	def __init__(self, probe, contacts, mapGraph, direction = True, path = []):
		Behavior.__init__(self, probe)

		self.mapGraph = mapGraph
		self.localNode = self.mapGraph.getCurrentNode()
		self.contacts = contacts
		self.frontAnchorFit = 0
		self.concertinaFit = 0

		
		self.direction = direction
		self.isInit = False

		self.holdT = HoldTransition(probe)
		self.holdSlideT = HoldSlideTransition(probe, direction)

		self.cPoints = []
		
		self.isJerking = False
		self.jerkingDone = False
		self.jerkAngle = 0
		self.prevJerkAngle = self.jerkAngle
		self.jerkErrors = []
		self.jerkJoint = 0
		self.jerkJoints = [0,1,2,3]
		self.nomJerk = 0.0
		self.nomJerks = [0.0,0.0,0.0,0.0]
		self.jerkAngles = [0.0,0.0,0.0,0.0]
		self.prevJerkAngles = [0.0,0.0,0.0,0.0]
		
		if direction:
			self.spliceJoint = 7
			#self.spliceJoint = 13
		else:
			self.spliceJoint = 31
			
		self.lastSpliceAngle = 0.0

		
		#self.computeCurve()

		self.mask = [0.0 for i in range(0,40)]
		self.count = 0
		
		self.currPeak = 0
		
		self.minAmp = 0.0
		self.maxAmp = 0.0
		self.ampInc = 0.04
			
		self.plotCount = 0
		
		self.newJoints = [None for i in range(0,self.probe.numSegs-1)]

		" flag true when the HoldTransition behavior has completed its behavior, select next transition "
		self.transDone = False

		" flag true when all the joints are stable, then step the HoldTransition behavior "
		self.refDone = False

		self.frontAnchoringState = True
		self.frontAnchoringDone = False
		self.frontExtending = True
		#self.frontExtending = False
		self.frontExtendDone = False
		

		self.globalCurveFit = 0
		#self.setPath(path)
				
	def getMask(self):
		return self.mask

	def setPath(self, path):
		print "setting path: ", path

		self.path = deepcopy(path)
		#self.computeCurve()

		self.pathCurve = VoronoiFit(self.path)
		
		self.globalCurveFit = GlobalCurveFit(self.probe, self.contacts, self.pathCurve)
		self.globalCurveSlide = GlobalCurveSlide(self.probe, self.contacts, self.pathCurve)
		
		#self.setDirection(not self.globalCurveFit.getPathDirection())
		result = self.globalCurveSlide.getPathDirection()
		
		print "received", result, "from globalCurveSlide"
		
		self.setDirection(result)
		
		#self.globalCurveFit.draw()
		self.globalCurveSlide.draw()
		
	def reverseDirection(self):
		self.path.reverse()
		
		self.setPath(self.path)
		
		self.computeCurve()
						
	def setDirection(self, isForward):

		self.direction = isForward
		
		print "setting direction =", self.direction

		self.mask = [0.0 for i in range(0,40)]
		self.count = 0

		self.frontAnchoringState = True
		self.frontAnchoringDone = False
		self.frontExtending = True
		self.frontExtendDone = False
					
		self.minAmp = 0.0
		self.maxAmp = 0.0			
		self.ampInc = 0.04
		
		self.currPeak = 0
		#self.spliceJoint = 7
		self.lastSpliceAngle = 0.0
						
		if self.direction:
			self.spliceJoint = 7
			#self.spliceJoint = 13
		else:
			self.spliceJoint = 31
		
		"reset everything "
		self.frontAnchorFit = 0
		self.frontCurve = 0
		self.adaptiveCurve = 0
		self.concertinaFit = 0
		
		#self.computeCurve()
		
	def computeCurve(self):

				
		self.frontCurve = AdaptiveAnchorCurve(4*pi, self.probe.segLength)				

		if self.frontAnchorFit == 0:
			self.frontAnchorFit = FrontAnchorFit(self.probe, self.direction, self.spliceJoint)
			self.frontAnchorFit.setCurve(self.frontCurve)
			self.frontAnchorFit.setSpliceJoint(self.spliceJoint)

		else:
			self.frontAnchorFit.setCurve(self.frontCurve)
			self.frontAnchorFit.setSpliceJoint(self.spliceJoint)


		self.adaptiveCurve = BackConcertinaCurve(4*pi)
		self.adaptiveCurve.setTailLength(2.0)
		
		if self.concertinaFit == 0:
			self.concertinaFit = FastLocalCurveFit(self.probe, self.direction, self.adaptiveCurve, 38)

			if self.direction:	
				self.concertinaFit.setStartNode(self.spliceJoint)
			else:
				self.concertinaFit.setEndNode(self.spliceJoint)
		else:
			self.concertinaFit.setCurve(self.adaptiveCurve)
			if self.direction:	
				self.concertinaFit.setStartNode(self.spliceJoint)
			else:
				self.concertinaFit.setEndNode(self.spliceJoint)


		self.frontAnchorFit.step()
		self.concertinaFit.step()
		
		resultJoints = self.spliceFitJoints()
		
		" straighten the part of frontAnchorFit that is not actuated "
		if self.frontAnchorFit.anterior:
			startNode = 0
			endNode = self.frontAnchorFit.lastJoint
			
			#self.globalCurveSlide.setBoundaries(startNode, 20)
			
			
		else:
			startNode = self.frontAnchorFit.lastJoint
			endNode = self.probe.numSegs-2
			#self.globalCurveSlide.setBoundaries(19, endNode)
			
		for i in range(startNode, endNode+1):
			print "changing joint", i, "from", resultJoints[i], "to 0.0"
			resultJoints[i] = 0.0
			

		
		#self.holdSlideT.setSpliceJoint(self.spliceJoint)

		self.globalCurveSlide.reset(resultJoints)
		#print "globalCurveSlide.step()"
		self.globalCurveSlide.step()
		resultJoints2 = self.globalCurveSlide.getJoints()
		
		#self.holdSlideT.reset(resultJoints, self.direction)
		self.holdT.reset(resultJoints2)
		self.refDone = False

	def spliceFitJoints(self, isStart = False):
		
		result = [None for i in range(self.probe.numSegs-1)]
		
		#print "spliceJoint = ", self.spliceJoint
		#print "frontAnchorJoints =", joints1
		#print "concertinaJoints =", joints2

		#if not self.concertinaFit.isJointSolid(self.spliceJoint):
		#	newJointVal = joints1[self.spliceJoint] + joints2[self.spliceJoint]
		#	self.lastSpliceAngle = newJointVal
		#else:
		#	newJointVal = self.lastSpliceAngle
				
		joints1 = self.frontAnchorFit.getJoints()
		joints2 = self.concertinaFit.getJoints()
		
		#print "joints1 =", joints1
		#print "joints2 =", joints2
		
		self.concertinaFit.spliceJoint = self.spliceJoint
		self.concertinaFit.spliceAngle = pi / 180.0 * joints1[self.spliceJoint]
				
		newJointVal = joints1[self.spliceJoint] + joints2[self.spliceJoint]
				
		tempJoints = copy(joints1)
		tempJoints[self.spliceJoint] = newJointVal
		
		jointList = [tempJoints,joints2]
		
		# merge the joints prioritizing the first to last
		for joints in jointList:
			for i in range(len(joints)):
				if result[i] == None and joints[i] != None:
					result[i] = joints[i]
		
		return result
	
	def getCenterPoints(self):
		
		#return copy(self.cmdSegPoints)
		return self.cPoints
		
	def computeCenterPoints(self):

		self.computeCmdSegPoints()
		
		#points = [[self.adaptiveCurve.ampSpacing/2.0 + i * self.adaptiveCurve.ampSpacing, -self.frontCurve.peakAmp[0]] for i in range(self.frontCurve.numPeaks)]
		
		points = [[0.0,0.0], [pi/self.frontCurve.initFreq, 0.0], [2*pi/self.frontCurve.initFreq, 0.0]]
		points.reverse()
		
		" 1. at each inflection point, find the closest segment and it's point on segment"
		" 2. cycle through all 40 segments "
		" 3. transform those points to current estimated positions "

		self.localNode = self.mapGraph.getCurrentNode()
		
		c_configs = []
		
		for p in points:
			minDist = 1e100
			minPoint = [0.0,0.0]
			jointIndex = 0
			segLen = 0.0
			
			for i in range(len(self.cmdSegPoints)-1):

				p0 = self.cmdSegPoints[i]
				p1 = self.cmdSegPoints[i+1]
				seg = [p0,p1]
				dist, cp = closestSegPoint(seg, p)

				if dist < minDist:
					minDist = dist
					minPoint = cp
					jointIndex = i-1
					segLen = sqrt((p0[0]-minPoint[0])**2 + (p0[1]-minPoint[1])**2)

			" return joint number and position on segment "
			c_configs.append([jointIndex,segLen])

		self.cPoints = []
		for c in c_configs:
			joint = c[0]
			segLen = c[1]
			
			pose = self.contacts.getAveragePose(joint)
			#pose = self.localNode.getJointPose(joint)
			point = [pose[0] + segLen*cos(pose[2]), pose[1] + segLen*sin(pose[2])]
			
			if point not in self.cPoints:
				self.cPoints.append(point)
		

	def computeCmdSegPoints(self):

		if self.direction:

			spliceJointAngle = self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0
						
			self.cmdSegPoints = []
			cmdOrigin = [0.0,0.0,-spliceJointAngle]
			self.cmdSegPoints.append(cmdOrigin)
			
			ind = range(0,self.spliceJoint+1)
			ind.reverse()
			
			for i in ind:		
				sampAngle = self.probe.getServoCmd(i)
				totalAngle = cmdOrigin[2] + sampAngle
				xTotal = cmdOrigin[0] - self.probe.segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] - self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.insert(0, cmdOrigin)

			cmdOrigin = [self.probe.segLength*cos(-spliceJointAngle),self.probe.segLength*sin(-spliceJointAngle),-spliceJointAngle]
			self.cmdSegPoints.append(cmdOrigin)

			for i in range(self.spliceJoint+1,self.probe.numSegs-1):
				
				sampAngle = self.probe.getServoCmd(i)

				totalAngle = cmdOrigin[2] - sampAngle
				xTotal = cmdOrigin[0] + self.probe.segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] + self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.append(cmdOrigin)
			
		else:
			
			spliceJointAngle = -self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0 
			
			self.cmdSegPoints = []
			cmdOrigin = [self.probe.segLength*cos(-spliceJointAngle),self.probe.segLength*sin(-spliceJointAngle),-spliceJointAngle + pi]

			self.cmdSegPoints.append(cmdOrigin)

			ind = range(0,self.spliceJoint)
			ind.reverse()
			
			for i in ind:
				sampAngle = self.probe.getServoCmd(i)
				totalAngle = cmdOrigin[2] + sampAngle
				xTotal = cmdOrigin[0] - self.probe.segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] - self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.insert(0,cmdOrigin)

			cmdOrigin = [0.0,0.0,-spliceJointAngle + pi]			
			self.cmdSegPoints.append(cmdOrigin)
			
			for i in range(self.spliceJoint,self.probe.numSegs-1):
				sampAngle = self.probe.getServoCmd(i)
				totalAngle = cmdOrigin[2] - sampAngle
				xTotal = cmdOrigin[0] + self.probe.segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] + self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.append(cmdOrigin)
					
	def drawFit(self):
		
		#peakJoints = self.frontAnchorFit.getPeakJoints(self.currPeak)

		pylab.clf()
		crvPoints = self.frontCurve.getPoints()		
		xP = []
		yP = []
		for p in crvPoints:
			px = p[0]*cos(pi) - p[1]*sin(pi)
			py = p[0]*sin(pi) + p[1]*cos(pi)
			
			xP.append(px)
			yP.append(py)
		pylab.plot(xP,yP,color='0.3')

		#pylab.clf()
		crvPoints = self.adaptiveCurve.getPoints()		
		xP = []
		yP = []
		for p in crvPoints:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP,color='0.3')
		
		self.computeCmdSegPoints()

		if self.direction:
			" the actual configuration of the snake "
			#actSegPoints = [[-self.probe.segLength,0.0,0.0]]

			#resultJoints = self.spliceFitJoints(self.frontAnchorFit.getJoints(), self.concertinaFit.getJoints())
			#spliceJointAngle = resultJoints[self.spliceJoint]
			#if not self.concertinaFit.isJointSolid(self.spliceJoint):
			#	spliceJointAngle = self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0
			#	self.lastDrawSpliceAngle = spliceJointAngle
			#else:
			#	spliceJointAngle = self.lastDrawSpliceAngle
	
			spliceJointAngle = self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0
			
			actSegPoints = []
			actOrigin = [0.0,0.0,-spliceJointAngle]

			actSegPoints.append(actOrigin)

			ind = range(0,self.spliceJoint+1)
			ind.reverse()
			
			for i in ind:
				sampAngle = self.probe.getServo(i)
				totalAngle = actOrigin[2] + sampAngle
				xTotal = actOrigin[0] - self.probe.segLength*cos(totalAngle)
				zTotal = actOrigin[1] - self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.insert(0,actOrigin)
	
			actOrigin = [self.probe.segLength*cos(-spliceJointAngle),self.probe.segLength*sin(-spliceJointAngle),-spliceJointAngle]
			actSegPoints.append(actOrigin)
			
			for i in range(self.spliceJoint+1,self.probe.numSegs-1):
				sampAngle = self.probe.getServo(i)
				totalAngle = actOrigin[2] - sampAngle
				xTotal = actOrigin[0] + self.probe.segLength*cos(totalAngle)
				zTotal = actOrigin[1] + self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.append(actOrigin)

		else:
			
			spliceJointAngle = -self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0 
			
			actSegPoints = []
			actOrigin = [self.probe.segLength*cos(-spliceJointAngle),self.probe.segLength*sin(-spliceJointAngle),-spliceJointAngle + pi]

			actSegPoints.insert(0, actOrigin)

			ind = range(0,self.spliceJoint)
			ind.reverse()
			
			for i in ind:
				sampAngle = self.probe.getServo(i)
				totalAngle = actOrigin[2] + sampAngle
				xTotal = actOrigin[0] - self.probe.segLength*cos(totalAngle)
				zTotal = actOrigin[1] - self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.insert(0,actOrigin)

			actOrigin = [0.0,0.0,-spliceJointAngle + pi]			
			actSegPoints.append(actOrigin)
			
			for i in range(self.spliceJoint,self.probe.numSegs-1):
				sampAngle = self.probe.getServo(i)
				totalAngle = actOrigin[2] - sampAngle
				xTotal = actOrigin[0] + self.probe.segLength*cos(totalAngle)
				zTotal = actOrigin[1] + self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.append(actOrigin)

			
			" the actual configuration of the snake "
			"""
			actSegPoints = [[-self.probe.segLength,0.0,0.0]]
			actOrigin = [0.0,0.0,pi]
			actSegPoints.insert(0,actOrigin)

			ind = range(0,self.probe.numSegs-1)
			ind.reverse()
				
			for i in ind:
				sampAngle = self.probe.getServo(i)
				totalAngle = actOrigin[2] + sampAngle
				xTotal = actOrigin[0] - self.probe.segLength*cos(totalAngle)
				zTotal = actOrigin[1] - self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.insert(0,actOrigin)
			
				
				#actSegPoints[i+2] = actOrigin		
			"""
		
		if self.direction:

			
			spliceJointAngle = self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0

			" local peak configurations "
			peakCurves = [[]]

			peakOrigin = [self.probe.segLength*cos(-spliceJointAngle),self.probe.segLength*sin(-spliceJointAngle),-spliceJointAngle]
			peakCurves[-1].append(copy(peakOrigin))


			for i in range(self.spliceJoint+1,self.probe.numSegs-1):
				sampAngle = self.probe.getServo(i)
				totalAngle = peakOrigin[2] - sampAngle
				xTotal = peakOrigin[0] + self.probe.segLength*cos(totalAngle)
				zTotal = peakOrigin[1] + self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				peakOrigin = pnt
				#peakCurves[-1].append(copy(peakOrigin))
				
				" if current joint is the minimum of a current peak, reset the origin "
				peakIndex = -1
				for index in range(len(self.concertinaFit.jointClasses)):
					if self.concertinaFit.jointClasses[index].count(i) > 0:
						peakIndex = index
						break
				
				if peakIndex >= 1 and min(self.concertinaFit.jointClasses[peakIndex]) == i:
					peakCurves.append([])
					peakOrigin = self.cmdSegPoints[i+1]
					peakCurves[-1].append(copy(peakOrigin))
	
					sampAngle = self.probe.getServo(i)
					totalAngle = peakOrigin[2] - sampAngle
					xTotal = peakOrigin[0] + self.probe.segLength*cos(totalAngle)
					zTotal = peakOrigin[1] + self.probe.segLength*sin(totalAngle)
					pnt = [xTotal, zTotal, totalAngle]
					peakOrigin = pnt
					peakCurves[-1].append(copy(peakOrigin))
				else:
					peakCurves[-1].append(copy(peakOrigin))
		else:

			spliceJointAngle = -self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0

			" local peak configurations "
			peakCurves = [[]]

			peakOrigin = [self.probe.segLength*cos(-spliceJointAngle),self.probe.segLength*sin(-spliceJointAngle),-spliceJointAngle + pi]
			peakCurves[-1].append(copy(peakOrigin))
			
			ind = range(0,self.spliceJoint)
			ind.reverse()

			for i in ind:
				sampAngle = self.probe.getServo(i)
				totalAngle = peakOrigin[2] + sampAngle
				xTotal = peakOrigin[0] - self.probe.segLength*cos(totalAngle)
				zTotal = peakOrigin[1] - self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				peakOrigin = pnt
				#peakCurves[-1].append(copy(peakOrigin))
				
				" if current joint is the minimum of a current peak, reset the origin "
				peakIndex = -1
				for index in range(len(self.concertinaFit.jointClasses)):
					if self.concertinaFit.jointClasses[index].count(i) > 0:
						peakIndex = index
						break
				
				if peakIndex >= 1 and max(self.concertinaFit.jointClasses[peakIndex]) == i:
					peakCurves.append([])
					peakOrigin = self.cmdSegPoints[i+1]
					peakCurves[-1].append(copy(peakOrigin))
	
					sampAngle = self.probe.getServo(i)
					totalAngle = peakOrigin[2] + sampAngle
					xTotal = peakOrigin[0] - self.probe.segLength*cos(totalAngle)
					zTotal = peakOrigin[1] - self.probe.segLength*sin(totalAngle)
					pnt = [xTotal, zTotal, totalAngle]
					peakOrigin = pnt
					peakCurves[-1].append(copy(peakOrigin))
				else:
					peakCurves[-1].append(copy(peakOrigin))
					
					
							

			"""
			spliceJointAngle = -self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0

			" local peak configurations "
			peakCurves = [[]]

			#peakOrigin = [self.probe.segLength*cos(-spliceJointAngle),self.probe.segLength*sin(-spliceJointAngle),-spliceJointAngle + pi]

			#peakOrigin = [self.probe.segLength*cos(-spliceJointAngle),self.probe.segLength*sin(-spliceJointAngle),-spliceJointAngle]
			
			peakOrigin = [0.0,0.0,-spliceJointAngle + pi]	
			peakCurves[-1].append(copy(peakOrigin))

			ind = range(0,self.spliceJoint)
			ind.reverse()
			
			#ind = range(0, self.spliceJoint+2)
			#ind.reverse()

			for i in ind:
				sampAngle = self.probe.getServo(i)
				totalAngle = peakOrigin[2] + sampAngle
				xTotal = peakOrigin[0] - self.probe.segLength*cos(totalAngle)
				zTotal = peakOrigin[1] - self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				peakOrigin = pnt
				
				" if current joint is the minimum of a current peak, reset the origin "
				peakIndex = -1
				for index in range(len(self.concertinaFit.jointClasses)):
					if self.concertinaFit.jointClasses[index].count(i) > 0:
						peakIndex = index
						break
				
				if peakIndex >= 1 and min(self.concertinaFit.jointClasses[peakIndex]) == i:
					peakCurves.append([])
					peakOrigin = self.cmdSegPoints[i+1]
					peakCurves[-1].append(copy(peakOrigin))
	
					sampAngle = self.probe.getServo(i)
					totalAngle = peakOrigin[2] + sampAngle
					xTotal = peakOrigin[0] - self.probe.segLength*cos(totalAngle)
					zTotal = peakOrigin[1] - self.probe.segLength*sin(totalAngle)
					pnt = [xTotal, zTotal, totalAngle]
					peakOrigin = pnt
					peakCurves[-1].append(copy(peakOrigin))
				else:
					peakCurves[-1].append(copy(peakOrigin))
			"""
			
			
		errors = []
		for i in range(0,self.probe.numSegs-2):
			errors.append(fabs(self.probe.getServo(i)-self.probe.getServoCmd(i)))

		if True:
			xP = []
			yP = []
			for p in actSegPoints:
				xP.append(p[0])
				yP.append(p[1])		
			pylab.plot(xP,yP,color='r')


		if True:
			xP = []
			yP = []
			for p in self.cmdSegPoints:
				xP.append(p[0])
				yP.append(p[1])		
			pylab.scatter(xP,yP,color='g')
		
		if True	:
			for curve in peakCurves:
				xP = []
				yP = []
				for p in curve:
					xP.append(p[0])
					yP.append(p[1])
					
					
				pylab.scatter(xP,yP, linewidth=1, color='b')
				pylab.plot(xP,yP,color='b')
		
		"""
		for i in range(1, len(actSegPoints)-2):
			xP = [actSegPoints[i][0]]
			yP = [actSegPoints[i][1]]
			val = errors[i-1]
			if val > 1.0:
				val = 1.0
				
			pylab.scatter(xP,yP,color=str(1.0-val))
		"""

		#print "plotting graph", self.plotCount

		pylab.xlim(-4,4)
		#pylab.xlim(-0.5,4)
		pylab.ylim(-3,3)
		pylab.savefig("fitPlot%04u.png" % self.plotCount)
		
		self.plotCount += 1
		
	def computeAnchorErrors(self):
		numPoses = len(self.errorPoses1)
		xDiff = 0.0
		yDiff = 0.0
		for p in self.errorPoses1:
			xDiff += p[0]
			yDiff += p[1]
			
		xDiff /= numPoses
		yDiff /= numPoses
		
		err1 = sqrt(xDiff**2 + yDiff**2)

		numPoses = len(self.errorPoses2)
		xDiff = 0.0
		yDiff = 0.0
		for p in self.errorPoses2:
			xDiff += p[0]
			yDiff += p[1]
			
		xDiff /= numPoses
		yDiff /= numPoses
		
		err2 = sqrt(xDiff**2 + yDiff**2)
		
		return err1, err2
			
	def doExtendFront(self):

		print "Extend Front"

		#self.frontExtendDone = False
		
		if self.globalCurveSlide.step():
			self.frontExtendDone = True
			self.frontExtending = False
			
		joints = self.globalCurveSlide.getJoints()

		self.holdT.reset(joints)
		
		#resultJoints = self.spliceFitJoints()
		#self.holdSlideT.reset(resultJoints, self.direction)
		

		#self.frontExtending = False
		#self.frontExtendDone = False
					
	def doFrontAnchor(self):
				
		" compute the maximum error of the joints in the current peak "
		anchorJoints = self.frontAnchorFit.getPeakJoints()
		anchorJoints.sort()

		errors = []
		for j in anchorJoints:
			errors.append(self.probe.getServo(j)-self.probe.getServoCmd(j))

		maxError = 0.0	
		for err in errors:
			if fabs(err) > maxError:
				maxError = fabs(err)
		
		" check for terminating criteria "
		currAmp = self.frontCurve.getPeakAmp()	
		
		" draw the state of the curve fitting "
		#self.drawFit()

		currAmp = self.frontCurve.getPeakAmp()
		nextVal = currAmp	
		
		if not self.isJerking:
			" perform amplitude operation here "

			print "amps before:", self.minAmp, self.maxAmp, self.ampInc, currAmp, nextVal
			
			" error threshold that determines we have made an anchor contact "
			if maxError > 0.3 and self.maxAmp != 0.0:
				
				#print "errors =", errors
			
				" if the difference between amplitudes is < 0.01 "
				if fabs(self.maxAmp-self.minAmp) < 0.01:
					
					" A: Section of snake successfully anchored "
					
					" set our new peak to solid so we don't recompute its values "
					if self.currPeak != 0:
						self.currPeak += 2
						self.minAmp = 0.0
						self.maxAmp = 0.0
						
						nextVal = 0.0
						self.ampInc = 0.04
	
					else:
						
						self.isJerking = True
						

				else:
				
					" B: Section of snake experiencing error, now narrow down a tight fit"
				
					currAmp = self.frontCurve.getPeakAmp()
					
					" FIXME:  will these comparisons always work because of floating point error? "							
					" FIXME:  if self.minAmp is 0 and currAmp is 0, will reduce to negative amplitude "
					" error with zero amplitude means that the straight segments are already under stress "
					if currAmp == self.minAmp:
						
						" C:  Overshoot.  Minimum is too high." 
						" Reduce amplitude minimum, set current to maximum"
						
						" FIXME:  sometimes the amplitude goes below zero. "
						self.minAmp -= self.ampInc
						self.maxAmp = currAmp
						nextVal = self.minAmp
						
					else:
						
						"D: Bring amplitude down to minimum, and lower the step size"
						
						" bring it down to the minimum again "
						nextVal = self.minAmp

						" cut down the step size "
						self.ampInc /= 2.0
				

			else:
				
				"E:  Not experience any error, increase the amplitude by step size "

				" error not maxed out, so set min to the current amplitude "
				self.minAmp = currAmp

				nextVal = self.minAmp + self.ampInc
				
				" maximum is the value we just set "
				self.maxAmp = nextVal
			
			print "amps after:", self.minAmp, self.maxAmp, self.ampInc, currAmp, nextVal
				
		" change the amplitude "
		#print "setting amp = " , nextVal		
		self.frontCurve.setPeakAmp(nextVal)
		
			
		if self.isJerking:
			
			anchorJoints = self.frontAnchorFit.getPeakJoints() 
			if len(anchorJoints) > 0:
				
				if self.direction:
					maxJoint = max(anchorJoints)
					self.jerkJoint = maxJoint + 1
					
					self.jerkJoints = []
					for k in range(0,4):
						jointNum = maxJoint + k + 1
						if jointNum > 38:
							break
						self.jerkJoints.append(jointNum)
						
				else:
					minJoint = min(anchorJoints)
					self.jerkJoint = minJoint - 1
					
					self.jerkJoints = []
					for k in range(0,4):
						jointNum = minJoint - k - 1
						if jointNum < 0:
							break
						self.jerkJoints.append(jointNum)

				#print "anchor joints =", anchorJoints
				#print "setting jerkJoints to", self.jerkJoints
				
				if self.jerkAngle == 0 and self.prevJerkAngle == 0:
					self.nomJerk = self.probe.getServo(self.jerkJoint) * 180.0 / pi
				
					for k in range(len(self.jerkJoints)):
						self.nomJerks[k] = self.probe.getServo(self.jerkJoints[k]) * 180.0 / pi
				
					#print "nominal = " , self.nomJerk, self.probe.getServo(self.jerkJoint)
					self.prevJerkAngle = self.jerkAngle
					self.jerkAngle = 60	
				
					for k in range(len(self.jerkJoints)):
						self.prevJerkAngles[k] = self.jerkAngles[k]
						self.jerkAngles[k] = 60
						
					self.originPoses = []
					for k in anchorJoints:
						self.originPoses.append(self.contacts.getAveragePose(k))

				elif self.jerkAngle == 60:
					
					#print "setting jerk joint to ", self.jerkAngle + self.nomJerk
					#print "jerk angle error = " , self.probe.getServo(self.jerkJoint)-self.probe.getServoCmd(self.jerkJoint)

					errs = []
					for k in range(len(self.jerkJoints)):
						errs.append(self.probe.getServo(self.jerkJoints[k])-self.probe.getServoCmd(self.jerkJoints[k]))
						self.prevJerkAngles[k] = self.jerkAngles[k]
					
					self.jerkErrors.append(self.probe.getServo(self.jerkJoint)-self.probe.getServoCmd(self.jerkJoint))
					self.prevJerkAngle = self.jerkAngle
					self.jerkAngle = -60

					for k in range(len(self.jerkJoints)):
						self.prevJerkAngles[k] = self.jerkAngles[k]
						self.jerkAngles[k] = -60

					self.errorPoses1 = []
					for k in range(len(anchorJoints)):
						tempPose = self.contacts.getAveragePose(anchorJoints[k])
						originPose = self.originPoses[k]
						self.errorPoses1.append([tempPose[0]-originPose[0],tempPose[1]-originPose[1],tempPose[2]-originPose[2]])
				
					#print self.errorPoses1
					
					self.originPoses = []
					for k in anchorJoints:
						self.originPoses.append(self.contacts.getAveragePose(k))
					
				elif self.jerkAngle == -60:
					#print "jerk angle error = " , self.probe.getServo(self.jerkJoint)-self.probe.getServoCmd(self.jerkJoint)
					self.jerkErrors.append(self.probe.getServo(self.jerkJoint)-self.probe.getServoCmd(self.jerkJoint))

					" error = (-0.09, 0.005) is good "
					" error = (-2.73, -1.99) is bad "
					self.prevJerkAngle = self.jerkAngle
					self.jerkAngle = 0
					
					for k in range(len(self.jerkJoints)):
						self.prevJerkAngles[k] = self.jerkAngles[k]
						self.jerkAngles[k] = 0

					self.errorPoses2 = []
					for k in range(len(anchorJoints)):
						tempPose = self.contacts.getAveragePose(anchorJoints[k])
						originPose = self.originPoses[k]
						self.errorPoses2.append([tempPose[0]-originPose[0],tempPose[1]-originPose[1],tempPose[2]-originPose[2]])
				
					#print self.errorPoses2
				
				else:
					self.prevJerkAngle = self.jerkAngle
					self.jerkingDone = True

					for k in range(len(self.jerkJoints)):
						self.prevJerkAngles[k] = self.jerkAngles[k]

				#print "setting jerk joint to ", self.jerkAngle + self.nomJerk
				#print "jerk angle error = " , self.probe.getServo(j)-self.probe.getServoCmd(j)

		if self.isJerking and self.jerkingDone:
			
			" if the anchor is not secure, lets try it again "
			err1, err2 = self.computeAnchorErrors()
			
			#print "anchor errors =", err1, err2
			
			if err1 > 0.1 or err2 > 0.1:

				" TODO " 
				" if the anchor failed, increment the splice joint higher "

			else:
				
				" Front anchoring is done! "
				" we're done with this behavior. go to concertina gait "
				self.frontAnchoringDone = True

			" reset everything to default values "	
			self.minAmp = 0.0
			self.maxAmp = 0.0

			nextVal = 0.0
			self.ampInc = 0.04
		
			self.isJerking = False
			self.jerkingDone = False
			
			self.jerkErrors = []


		" reset all torques of the joints to maximum "
		for i in range(self.probe.numSegs-1):
			self.probe.setJointTorque(i, self.probe.maxTorque)

		" weaken head joints "
		" 30*2.5 is maximum "
		anchorJoints = self.frontAnchorFit.getPeakJoints()

		#print "anchorJoints =", len(anchorJoints), "amp =", self.frontCurve.getPeakAmp()
		#if len(anchorJoints) > 0 and self.frontCurve.getPeakAmp() > 0 and self.isJerking:

		" Lead joints of the snake are weakened when not part of front-anchoring "		
		if len(anchorJoints) > 0 and self.frontCurve.getPeakAmp() > 0:
			
			if not self.direction:
				maxJoint = max(anchorJoints)
				#print "weakening joints", range(maxJoint+1, self.probe.numSegs-2)
				for i in range(maxJoint+1, self.probe.numSegs-1):
					if i <= self.probe.numSegs-2:
						self.probe.setJointTorque(i, 3.0)
			else:
				minJoint = min(anchorJoints)
				#print "weakening joints", range(0,minJoint-1)
				for i in range(0, minJoint-1):					
					if i <= self.probe.numSegs-2:
						self.probe.setJointTorque(i, 3.0)
	
		" execute the local curve fitting "
		self.frontAnchorFit.step()

		" compute the bounds for behavior for global curve fitting "
		if self.frontAnchorFit.anterior:
			startNode = 0
			endNode = self.frontAnchorFit.lastJoint
			self.globalCurveFit.setBoundaries(startNode, endNode)
		else:
			startNode = self.frontAnchorFit.lastJoint
			endNode = self.probe.numSegs-2	
			self.globalCurveFit.setBoundaries(startNode, endNode)
			
		" Set the joints, with joints2 as priority over joints1 "
		#self.mergeJoints([joints2, joints1])

		resultJoints = self.spliceFitJoints()
		

		self.holdT.reset(resultJoints)
		
		if self.isJerking:
			#print "setting jerk joints", self.jerkJoints, "to", self.nomJerk + self.jerkAngle
			#self.holdT.positions[self.jerkJoint] = self.nomJerk + self.jerkAngle
			for k in range(len(self.jerkJoints)):
				self.holdT.positions[self.jerkJoints[k]] = self.nomJerks[k] + self.jerkAngles[k]
			
	def doBackConcertina(self):
		
		self.transDone = False
		
		" compute the maximum error of the joints in the current peak "
		peakJoints = self.concertinaFit.getPeakJoints(self.currPeak) + self.concertinaFit.getPeakJoints(self.currPeak+1)
		peakJoints.sort()
		
		errors = []
		for j in peakJoints:
			errors.append(self.probe.getServo(j)-self.probe.getServoCmd(j))

		maxError = 0.0	
		for err in errors:
			if fabs(err) > maxError:
				maxError = fabs(err)
		
		#print "maxError =", maxError
		
		" draw the state of the curve fitting "
		#self.drawFit()

		currAmp = self.adaptiveCurve.getPeakAmp(self.currPeak)
		nextVal = currAmp

		" perform amplitude operation here "
		
		" error threshold that determines we have made an anchor contact "
		if maxError > 0.3 and self.maxAmp != 0.0:
			
			" if the difference between amplitudes is < 0.01 "
			if fabs(self.maxAmp-self.minAmp) < 0.01:
				
				" A: Section of snake successfully anchored "
				
				" set our new peak to solid so we don't recompute its values "
				self.concertinaFit.setSolid(self.currPeak)
				self.concertinaFit.setSolid(self.currPeak+1)
				self.currPeak += 2
				self.minAmp = 0.0
				self.maxAmp = 0.0
				
				nextVal = 0.0
				self.ampInc = 0.04
				
			else:
			
				" B: Section of snake experiencing error, now narrow down a tight fit"
			
				currAmp = self.adaptiveCurve.getPeakAmp(self.currPeak)
				
				if currAmp == self.minAmp:
					
					" C:  Overshoot.  Minimum is too high." 
					" Reduce amplitude minimum, set current to maximum"
					
					self.minAmp -= self.ampInc
					self.maxAmp = currAmp
					nextVal = self.minAmp

				else:
					
					"D: Bring amplitude down to minimum, and lower the step size"
					
					" bring it down to the minimum again "
					nextVal = self.minAmp

					" cut down the step size "
					self.ampInc /= 2.0
			

		else:
			
			"E:  Not experience any error, increase the amplitude by step size "

			" error not maxed out, so set min to the current amplitude "
			self.minAmp = currAmp

			nextVal = self.minAmp + self.ampInc
			
			" maximum is the value we just set "
			self.maxAmp = nextVal


		" reset all torques of the joints to maximum "
		for i in range(self.probe.numSegs-1):
			self.probe.setJointTorque(i, self.probe.maxTorque)

		" stretch tail to remove closed-chain interference "
		self.adaptiveCurve.setTailLength(5.0)
		self.adaptiveCurve.setPeakAmp(self.currPeak,nextVal)
		self.adaptiveCurve.setPeakAmp(self.currPeak+1,nextVal)

		" execute the local curve fitting "
		self.concertinaFit.step()

		" weaken the feeder joints "
		" 30*2.5 is maximum "
		peakJoints = self.concertinaFit.getPeakJoints(self.currPeak+1)
		if len(peakJoints) > 0:
			
			if self.direction:
				maxJoint = max(peakJoints)
				#print "weakening", maxJoint+1, "to", self.probe.numSegs-3
				#for i in range(maxJoint+1, self.probe.numSegs-2):
				for i in range(maxJoint+2, self.probe.numSegs-1):
					if i <= self.probe.numSegs-2:
						self.probe.setJointTorque(i, 3.0)
			else:
				minJoint = min(peakJoints)
				#print "weakening", 0, "to", minJoint-1
				#for i in range(0, minJoint):					
				for i in range(0, minJoint-1):					
					if i <= self.probe.numSegs-2:
						self.probe.setJointTorque(i, 3.0)
							

		resultJoints = self.spliceFitJoints()
		self.holdT.reset(resultJoints)

		peakJoints = self.concertinaFit.getPeakJoints(self.currPeak) + self.concertinaFit.getPeakJoints(self.currPeak+1)
		peakJoints.sort()
		
		#print "amp =", nextVal, "min =", self.minAmp, "max =", self.maxAmp
		#print "peak =", self.currPeak
		#print "peakJoints =", peakJoints
		

	def step(self):
		
		self.count += 1
		isDone = False
		
		if self.transDone:
			self.transDone = False


			#if self.frontExtending:
				#self.frontExtendDone = True
				#self.frontExtending = False
				#resultJoints = self.holdSlideT.getJoints()
				#self.holdT.reset(resultJoints)
				
		
			# termination cases
			if not self.frontAnchoringState:
				
				peakJoints = self.concertinaFit.getPeakJoints(self.currPeak) + self.concertinaFit.getPeakJoints(self.currPeak+1)
				peakJoints.sort()
	
				currAmp = self.adaptiveCurve.getPeakAmp(self.currPeak)
	
				if len(peakJoints) == 0 and self.currPeak != 0:
					isDone = True
					print "terminating caseA"
					
				elif self.direction and 38 in peakJoints and currAmp > 0.5:
					isDone = True
					print "terminating caseB"
					
				elif not self.direction and 0 in peakJoints and currAmp > 0.5:
					isDone = True
					print "terminating caseC"	

			if isDone:
					
				" prevent the back anchor from slithering around or buckling "
				if self.direction:
					for i in range(self.probe.numSegs-8, self.probe.numSegs-2):				
						self.holdT.positions[i] = 180.0 / pi * self.probe.getServo(i)
				else:
					for i in range(0, 6):					
						self.holdT.positions[i] = 180.0 / pi * self.probe.getServo(i)
	
				self.holdT.step()
				
			else:
							
				if self.frontAnchoringState:
					
					if self.frontExtending:
						
						print "Front Extend Step"
						self.doExtendFront()
						#self.frontExtending = False
						#self.frontExtendDone = False
		
						if self.frontExtendDone:
							self.frontExtending = False
						
					else:
						print "Front Anchor Step"
						self.doFrontAnchor()
		
						if self.frontAnchoringDone:
							self.frontAnchoringDone = False
							self.frontAnchoringState = False
							self.frontExtended = False
	
				else:
					peakJoints = self.concertinaFit.getPeakJoints(self.currPeak)	
					
					print "Concertina Step"
					self.doBackConcertina()


				solids = []
				for i in range(self.probe.numSegs-1):
					solids.append(self.concertinaFit.isJointSolid(i))
				#print "solids =", solids

				torques = []
				for i in range(0,self.probe.numSegs-2):
					torques.append(self.probe.getJointTorque(i))
				
				#print "joints =", self.spliceFitJoints()

				#print torques

				"""
				if self.frontExtending:
					self.refDone = False
					print "caseA"
					self.transDone = self.holdSlideT.step()
				
				else:
				"""	
				
				self.refDone = False
				self.transDone = self.holdT.step()
				
		else:
			
			" make no movements until all reference nodes are activated "
			if self.refDone:
				"""
				if self.frontExtending:
					self.transDone = self.holdSlideT.step()
				else:
					self.transDone = self.holdT.step()					
				"""
				
				self.transDone = self.holdT.step()					
				
				" execute the global curve fitting "
				#if self.frontAnchoringState and not self.frontExtending:
				#	self.globalCurveFit.step()
		
		anchorJoints = self.frontAnchorFit.getJoints()
		
		#print self.refDone, self.transDone, self.frontAnchoringState, self.frontExtending, self.frontExtendDone
		
		"""
		if self.frontExtending:
			joints = self.holdSlideT.getJoints()
		else:
			joints = self.holdT.getJoints()
		"""

		joints = self.holdT.getJoints()
		
		if self.frontAnchoringState and self.frontExtendDone:
			" collect joint settings for both behaviors "
			#joints2 = self.globalCurveFit.getJoints()
			
			joints2 = [0.0 for i in range(40)]
			
			
			if self.frontAnchorFit.anterior:			
				for p in range(0,self.frontAnchorFit.lastJoint+1):
					joints[p] = joints2[p]
			else:
				for p in range(self.frontAnchorFit.lastJoint, self.probe.numSegs-1):
					joints[p] = joints2[p]

		" Set the joints "
		if self.refDone:

			self.mergeJoints([joints])
			
			" update mask for ContactReference "
			for i in range(self.probe.numSegs-1):
				#print self.frontAnchoringState, anchorJoints[i], self.concertinaFit.isJointSolid(i), joints[i]
				if (not self.frontAnchoringState and anchorJoints[i] != None) or self.concertinaFit.isJointSolid(i) or joints[i] == None:
					self.mask[i] = 1.0
				else:
					self.mask[i] = 0.0

		else:
			" transition the masks, and wait for them to be stable, then change the joints "

			allActive = True
			" update mask for ContactReference "
			for i in range(self.probe.numSegs-1):
				if (not self.frontAnchoringState and anchorJoints[i] != None) or self.concertinaFit.isJointSolid(i) or joints[i] == None:
					self.mask[i] = 1.0

					" check if all reference points have been created "
					if not self.contacts.activeRef[i]:
						allActive = False
					
			if allActive:
				self.refDone = True
		
		#print self.mask
		
		if self.frontCurve.isOutOfSegments(self.frontAnchorFit.lastPosition):

			if self.direction:
				self.spliceJoint += 2
			else:
				self.spliceJoint -= 2
				
			#print "out of segments, switching splice joint to", self.spliceJoint
			self.mask = [0.0 for i in range(0,40)]
			self.count = 0
			
			self.minAmp = 0.0
			self.maxAmp = 0.0			
			self.ampInc = 0.04		
			self.currPeak = 0
			
			self.lastSpliceAngle = 0.0
			self.frontAnchorFit = 0
			self.frontCurve = 0
			self.adaptiveCurve = 0
			self.concertinaFit = 0

			self.computeCurve()
			
		if isDone:
			
			self.computeCenterPoints()
			self.mask = [0.0 for i in range(0,40)]
			self.count = 0

			self.frontAnchoringState = True
			self.frontAnchoringDone = False
			self.frontExtending = True
			self.frontExtendDone = False
			
			self.minAmp = 0.0
			self.maxAmp = 0.0			
			self.ampInc = 0.04
						
			self.currPeak = 0
			#self.spliceJoint = 7
			self.lastSpliceAngle = 0.0
							
			if self.direction:
				self.spliceJoint = 7
				#self.spliceJoint = 13
			else:
				self.spliceJoint = 31
			
			"reset everything "
			self.frontAnchorFit = 0
			self.frontCurve = 0
			self.adaptiveCurve = 0
			self.concertinaFit = 0

			self.computeCurve()
			

		return isDone
	

	def reset(self):

		#self.computeCenterPoints()
		#self.frontCurve.saveAmps()
		
		"reset everything "
		
		self.frontAnchorFit = 0
		self.frontCurve = 0
		self.computeCurve()
		
		self.mask = [0.0 for i in range(0,40)]
		self.count = 0
	
		self.currPeak = 0
	
	
	
	
