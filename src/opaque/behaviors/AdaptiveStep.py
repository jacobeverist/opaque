
from functions import *
from copy import copy, deepcopy
import pylab
from math import pi, sin, cos, sqrt, fabs, acos
from Behavior import Behavior
from FrontAnchorFit import FrontAnchorFit
from FastLocalCurveFit import FastLocalCurveFit
from HoldTransition import HoldTransition
from HoldSlideTransition import HoldSlideTransition
from AdaptiveAnchorCurve import AdaptiveAnchorCurve
from BackConcertinaCurve import BackConcertinaCurve

"""
TODO

1.  Check condition in FrontAnchorFit when curve goes beyond joint 0
2.  If joint 0 passed and anchor is loose, start over with higher splice joint
3.  Changing splice joint will shorten tail of BackConcertinaCurve
4.  Add the front anchoring state and then the concertina gait state to behavior control

"""

class AdaptiveStep(Behavior):

	def __init__(self, robotParam, probeState, contacts, mapGraph, direction = True):
		Behavior.__init__(self, robotParam)

		self.probeState = probeState
		self.mapGraph = mapGraph
		self.localNode = self.mapGraph.getCurrentNode()
		self.contacts = contacts
		self.frontAnchorFit = 0
		self.concertinaFit = 0
		
		self.segLength = self.robotParam['segLength']
		self.numJoints = self.robotParam['numJoints'] 
		self.maxTorque = self.robotParam['maxTorque']
		
		self.direction = direction
		self.isInit = False
		
		self.lastAttempt = False
		self.holdT = HoldTransition(robotParam)
		self.holdSlideT = HoldSlideTransition(robotParam, direction)

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
		else:
			self.spliceJoint = 31
			
		self.topJoint = 19

		self.lastSpliceAngle = 0.0

		self.mask = [0.0 for i in range(0,40)]
		self.count = 0
		
		self.currPeak = 0
		
		self.minAmp = 0.0
		self.maxAmp = 0.0
		self.ampInc = 0.04
			
		self.plotCount = 0

		self.torques = [None for i in range(0,self.numJoints)]
		self.newJoints = [None for i in range(0,self.numJoints)]

		" flag true when the HoldTransition behavior has completed its behavior, select next transition "
		self.transDone = False

		" flag true when all the joints are stable, then step the HoldTransition behavior "
		self.refDone = False

		self.frontAnchoringState = True
		self.frontAnchoringDone = False
		self.frontExtending = True
		self.frontExtendDone = False

		self.computeCurve()
		
	def setTopJoint(self, topJoint):
		self.topJoint = topJoint
				
	def getMask(self):
		#print "return mask:", self.mask
		return self.mask
		
	def setDirection(self, isForward):

		self.direction = isForward

		if self.direction:
			self.spliceJoint = 7
		else:
			self.spliceJoint = 31
		
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
		else:
			self.spliceJoint = 31
		
		"reset everything "
		self.frontAnchorFit = 0
		self.frontCurve = 0
		self.adaptiveCurve = 0
		self.concertinaFit = 0
		
		self.computeCurve()
		
	def computeCurve(self):
		
		self.frontCurve = AdaptiveAnchorCurve(4*pi, self.segLength)				

		if self.frontAnchorFit == 0:
			self.frontAnchorFit = FrontAnchorFit(self.robotParam, self.direction, self.spliceJoint)
			self.frontAnchorFit.setCurve(self.frontCurve)
			self.frontAnchorFit.setSpliceJoint(self.spliceJoint)

		else:
			self.frontAnchorFit.setCurve(self.frontCurve)
			self.frontAnchorFit.setSpliceJoint(self.spliceJoint)


		self.adaptiveCurve = BackConcertinaCurve(4*pi)
		#self.adaptiveCurve.setTailLength(2.0)

		if self.concertinaFit == 0:
			self.concertinaFit = FastLocalCurveFit(self.robotParam, self.direction, self.adaptiveCurve, 38)

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


		#self.adaptiveCurve.setPeakAmp(0,0.2)
		#self.adaptiveCurve.setPeakAmp(1,0.2)
		#self.adaptiveCurve.setPeakAmp(4,0.33)
		self.frontAnchorFit.step(self.probeState)
		self.concertinaFit.step(self.probeState)
		
		resultJoints = self.spliceFitJoints()
		
		self.holdSlideT.setSpliceJoint(self.spliceJoint)		
		self.holdSlideT.reset(self.probeState, resultJoints, self.direction)
		
		self.holdT.reset(self.probeState, resultJoints)
		self.refDone = False

		self.computeMaskAndOutput()
		
		#self.adaptiveCurve.setTailLength(2.0)

	def spliceFitJoints(self, isStart = False):
		
		result = [None for i in range(self.numJoints)]
		
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
	
		if self.direction:
			joints1[0] = 0.0
			joints1[1] = 0.0
			joints1[2] = 0.0
		else:
			joints1[-1] = 0.0
			joints1[-2] = 0.0
			joints1[-3] = 0.0
			
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

		cmdJoints = self.probeState['cmdJoints']
		
		if self.direction:
			
			" the configuration of the snake as we've commanded it "
			#self.cmdSegPoints = [[-self.segLength,0.0,0.0]]
	
			#if not self.concertinaFit.isJointSolid(self.spliceJoint):
			#	spliceJointAngle = self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0
			#	self.lastDrawSpliceAngle = spliceJointAngle
			#else:
			#	spliceJointAngle = self.lastDrawSpliceAngle

			spliceJointAngle = self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0
			
			#spliceJointAngle = self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0
			#spliceJoint = 0.0
			
			#print "joint =", spliceJointAngle
			
			self.cmdSegPoints = []
			cmdOrigin = [0.0,0.0,-spliceJointAngle]
			self.cmdSegPoints.append(cmdOrigin)
			
			ind = range(0,self.spliceJoint+1)
			ind.reverse()
			
			for i in ind:		
				sampAngle = cmdJoints[i]
				totalAngle = cmdOrigin[2] + sampAngle
				xTotal = cmdOrigin[0] - self.segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] - self.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.insert(0, cmdOrigin)

			cmdOrigin = [self.segLength*cos(-spliceJointAngle),self.segLength*sin(-spliceJointAngle),-spliceJointAngle]
			self.cmdSegPoints.append(cmdOrigin)

			for i in range(self.spliceJoint+1,self.numJoints):
				
				sampAngle = cmdJoints[i]

				totalAngle = cmdOrigin[2] - sampAngle
				xTotal = cmdOrigin[0] + self.segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] + self.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.append(cmdOrigin)
			
			
		else:
			
			spliceJointAngle = -self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0 
			
			self.cmdSegPoints = []
			cmdOrigin = [self.segLength*cos(-spliceJointAngle),self.segLength*sin(-spliceJointAngle),-spliceJointAngle + pi]

			self.cmdSegPoints.append(cmdOrigin)

			ind = range(0,self.spliceJoint)
			ind.reverse()
			
			for i in ind:
				sampAngle = cmdJoints[i]
				totalAngle = cmdOrigin[2] + sampAngle
				xTotal = cmdOrigin[0] - self.segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] - self.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.insert(0,cmdOrigin)

			cmdOrigin = [0.0,0.0,-spliceJointAngle + pi]			
			self.cmdSegPoints.append(cmdOrigin)
			
			for i in range(self.spliceJoint,self.numJoints):
				sampAngle = cmdJoints[i]
				totalAngle = cmdOrigin[2] - sampAngle
				xTotal = cmdOrigin[0] + self.segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] + self.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.append(cmdOrigin)
			
	def drawFit(self):

		cmdJoints = self.probeState['cmdJoints']
		
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

		stateJoints = self.probeState['joints']

		if self.direction:
			" the actual configuration of the snake "
			#actSegPoints = [[-self.segLength,0.0,0.0]]

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
				sampAngle = stateJoints[i]
				totalAngle = actOrigin[2] + sampAngle
				xTotal = actOrigin[0] - self.segLength*cos(totalAngle)
				zTotal = actOrigin[1] - self.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.insert(0,actOrigin)
	
			actOrigin = [self.segLength*cos(-spliceJointAngle),self.segLength*sin(-spliceJointAngle),-spliceJointAngle]
			actSegPoints.append(actOrigin)
			
			for i in range(self.spliceJoint+1,):
				sampAngle = stateJoints[i]
				totalAngle = actOrigin[2] - sampAngle
				xTotal = actOrigin[0] + self.segLength*cos(totalAngle)
				zTotal = actOrigin[1] + self.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.append(actOrigin)

		else:
			
			spliceJointAngle = -self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0 
			
			actSegPoints = []
			actOrigin = [self.segLength*cos(-spliceJointAngle),self.segLength*sin(-spliceJointAngle),-spliceJointAngle + pi]

			actSegPoints.insert(0, actOrigin)

			ind = range(0,self.spliceJoint)
			ind.reverse()
			
			for i in ind:
				sampAngle = stateJoints[i]
				totalAngle = actOrigin[2] + sampAngle
				xTotal = actOrigin[0] - self.segLength*cos(totalAngle)
				zTotal = actOrigin[1] - self.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.insert(0,actOrigin)

			actOrigin = [0.0,0.0,-spliceJointAngle + pi]			
			actSegPoints.append(actOrigin)
			
			for i in range(self.spliceJoint,self.numJoints):
				sampAngle = stateJoints[i]
				totalAngle = actOrigin[2] - sampAngle
				xTotal = actOrigin[0] + self.segLength*cos(totalAngle)
				zTotal = actOrigin[1] + self.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.append(actOrigin)

		if self.direction:

			
			spliceJointAngle = self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0

			" local peak configurations "
			peakCurves = [[]]

			peakOrigin = [self.segLength*cos(-spliceJointAngle),self.segLength*sin(-spliceJointAngle),-spliceJointAngle]
			peakCurves[-1].append(copy(peakOrigin))


			for i in range(self.spliceJoint+1,self.numJoints):
				sampAngle = stateJoints[i]
				totalAngle = peakOrigin[2] - sampAngle
				xTotal = peakOrigin[0] + self.segLength*cos(totalAngle)
				zTotal = peakOrigin[1] + self.segLength*sin(totalAngle)
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
	
					sampAngle = stateJoints[i]
					totalAngle = peakOrigin[2] - sampAngle
					xTotal = peakOrigin[0] + self.segLength*cos(totalAngle)
					zTotal = peakOrigin[1] + self.segLength*sin(totalAngle)
					pnt = [xTotal, zTotal, totalAngle]
					peakOrigin = pnt
					peakCurves[-1].append(copy(peakOrigin))
				else:
					peakCurves[-1].append(copy(peakOrigin))
		else:

			spliceJointAngle = -self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0

			" local peak configurations "
			peakCurves = [[]]

			peakOrigin = [self.segLength*cos(-spliceJointAngle),self.segLength*sin(-spliceJointAngle),-spliceJointAngle + pi]
			peakCurves[-1].append(copy(peakOrigin))
			
			ind = range(0,self.spliceJoint)
			ind.reverse()

			for i in ind:
				sampAngle = stateJoints[i]
				totalAngle = peakOrigin[2] + sampAngle
				xTotal = peakOrigin[0] - self.segLength*cos(totalAngle)
				zTotal = peakOrigin[1] - self.segLength*sin(totalAngle)
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
	
					sampAngle = stateJoints[i]
					totalAngle = peakOrigin[2] + sampAngle
					xTotal = peakOrigin[0] - self.segLength*cos(totalAngle)
					zTotal = peakOrigin[1] - self.segLength*sin(totalAngle)
					pnt = [xTotal, zTotal, totalAngle]
					peakOrigin = pnt
					peakCurves[-1].append(copy(peakOrigin))
				else:
					peakCurves[-1].append(copy(peakOrigin))
					
				
			
		errors = []
		for i in range(0,self.numJoints):
			errors.append(fabs(stateJoints[i]-cmdJoints[i]))

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

		print "plotting graph", self.plotCount

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

		self.transDone = self.holdSlideT.step(self.probeState)
		
		if self.transDone:
			self.frontExtendDone = True
			self.frontExtending = False

		print "doExtendFront"
			
	def doFrontAnchor(self):

		stateJoints = self.probeState['joints']
		cmdJoints = self.probeState['cmdJoints']

				
		" compute the maximum error of the joints in the current peak "
		anchorJoints = self.frontAnchorFit.getPeakJoints()
		anchorJoints.sort()

		errors = self.probeState['errors']

		maxError = 0.0	
		for j in anchorJoints:
			err = errors[j]
			if fabs(err) > maxError:
				maxError = fabs(err)
		
		
		" draw the state of the curve fitting "
		#self.drawFit()

		" check for terminating criteria "
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
						
						" do not allow the amplitude to go negative "
						if self.minAmp < 0.0:
							self.minAmp = 0.0
												
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
					self.nomJerk = stateJoints[self.jerkJoint] * 180.0 / pi
				
					for k in range(len(self.jerkJoints)):
						self.nomJerks[k] = stateJoints[self.jerkJoints[k]] * 180.0 / pi
				
					#print "nominal = " , self.nomJerk, stateJoints[self.jerkJoint)
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
					#print "jerk angle error = " , stateJoints[self.jerkJoint]-cmdJoints[self.jerkJoint]

					errs = []
					for k in range(len(self.jerkJoints)):
						errs.append(stateJoints[self.jerkJoints[k]]-cmdJoints[self.jerkJoints[k]])
						self.prevJerkAngles[k] = self.jerkAngles[k]
					
					self.jerkErrors.append(stateJoints[self.jerkJoint]-cmdJoints[self.jerkJoint])
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
					#print "jerk angle error = " , stateJoints[self.jerkJoint]-cmdJoints[self.jerkJoint]
					self.jerkErrors.append(stateJoints[self.jerkJoint]-cmdJoints[self.jerkJoint])

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
				#print "jerk angle error = " , stateJoints[j]-cmdJoints[j]

		if self.isJerking and self.jerkingDone:
			
			" if the anchor is not secure, lets try it again "
			err1, err2 = self.computeAnchorErrors()
			
			print "errorPoses1:", self.errorPoses1
			print "errorPoses2:", self.errorPoses2
			print "anchor errors =", err1, err2
			
			if err1 > 0.1 or err2 > 0.1:

				" TODO " 
				" if the anchor failed, increment the splice joint higher "

				" if we've exhausted too many joints, lets just quit with what we have "
				#if self.direction and self.spliceJoint >= 15:
				#	self.frontAnchoringDone = True
				#elif not self.direction and self.spliceJoint <= 23:
				#	self.frontAnchoringDone = True
				
				if self.lastAttempt:
					self.frontAnchoringDone = True
					print "ERROR: Quit with error on jerk joints in Front Anchoring"
					
				else:
					if self.direction:
						self.spliceJoint += 2
						
						if self.spliceJoint >= self.topJoint:
							self.spliceJoint = self.topJoint
							self.lastAttempt = True

					else:
						self.spliceJoint -= 2
						if self.spliceJoint <= self.topJoint:
							self.spliceJoint = self.topJoint
							self.lastAttempt = True
						
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
		for i in range(self.numJoints):
			self.torques[i] =  self.maxTorque

		" weaken head joints "
		" 30*2.5 is maximum "
		anchorJoints = self.frontAnchorFit.getPeakJoints()

		#print "anchorJoints =", len(anchorJoints), "amp =", self.frontCurve.getPeakAmp()
		#if len(anchorJoints) > 0 and self.frontCurve.getPeakAmp() > 0 and self.isJerking:

		" Lead joints of the snake are weakened when not part of front-anchoring "
		if len(anchorJoints) > 0 and self.frontCurve.getPeakAmp() > 0:
			
			if not self.direction:
				maxJoint = max(anchorJoints)
				for i in range(maxJoint+1, self.numJoints):
					if i <= self.numJoints-1:
						self.torques[i] = 3.0
			else:
				minJoint = min(anchorJoints)
				#print "weakening joints", range(0,minJoint-1)
				for i in range(0, minJoint-1):					
					if i <= self.numJoints-1:
						self.torques[i] = 3.0

		" execute the local curve fitting "
		self.frontAnchorFit.step(self.probeState)

		resultJoints = self.spliceFitJoints()
		self.holdT.reset(self.probeState, resultJoints)
		
		if self.isJerking:
			#print "setting jerk joints", self.jerkJoints, "to", self.nomJerk + self.jerkAngle
			#self.holdT.positions[self.jerkJoint] = self.nomJerk + self.jerkAngle
			for k in range(len(self.jerkJoints)):
				self.holdT.positions[self.jerkJoints[k]] = self.nomJerks[k] + self.jerkAngles[k]
		
			
	def doBackConcertina(self):
		
		stateJoints = self.probeState['joints']
		cmdJoints = self.probeState['cmdJoints']

		self.transDone = False
		
		" compute the maximum error of the joints in the current peak "
		peakJoints = self.concertinaFit.getPeakJoints(self.currPeak) + self.concertinaFit.getPeakJoints(self.currPeak+1)
		peakJoints.sort()
		
		errors = self.probeState['errors']
				
		maxError = 0.0	
		for j in peakJoints:
			err = errors[j]
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
					
					" do not allow the amplitude to go negative "
					if self.minAmp < 0.0:
						self.minAmp = 0.0
						
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
		for i in range(self.numJoints):
			self.torques[i] = self.maxTorque

		" stretch tail to remove closed-chain interference "
		self.adaptiveCurve.setTailLength(5.0)
		self.adaptiveCurve.setPeakAmp(self.currPeak,nextVal)
		self.adaptiveCurve.setPeakAmp(self.currPeak+1,nextVal)

		" execute the local curve fitting "
		self.concertinaFit.step(self.probeState)

		" weaken the feeder joints "
		" 30*2.5 is maximum "
		peakJoints = self.concertinaFit.getPeakJoints(self.currPeak+1)
		if len(peakJoints) > 0:
			
			if self.direction:
				maxJoint = max(peakJoints)
				for i in range(maxJoint+2, self.numJoints):
					if i <= self.numJoints-1:
						self.torques[i] = 3.0
			else:
				minJoint = min(peakJoints)
				#print "weakening", 0, "to", minJoint-1
				#for i in range(0, minJoint):					
				for i in range(0, minJoint-1):					
					if i <= self.numJoints-1:
						self.torques[i] = 3.0
							

		resultJoints = self.spliceFitJoints()
		self.holdT.reset(self.probeState, resultJoints)

		peakJoints = self.concertinaFit.getPeakJoints(self.currPeak) + self.concertinaFit.getPeakJoints(self.currPeak+1)
		peakJoints.sort()
		
		#print "amp =", nextVal, "min =", self.minAmp, "max =", self.maxAmp
		#print "peak =", self.currPeak
		#print "peakJoints =", peakJoints
		

	def step(self, probeState):
		Behavior.step(self, probeState)
		
		self.count += 1
		isDone = False

		stateJoints = self.probeState['joints']
		torques = self.probeState['torques']

		self.torques = []
		for i in range(0,self.numJoints):
			self.torques.append(torques[i])
			
		#print self.count, self.transDone, self.refDone, self.frontAnchoringState, self.frontExtending, self.frontExtendDone
		
		if self.transDone:
			self.transDone = False

			if self.frontExtending:
				self.frontExtendDone = True
				self.frontExtending = False
				resultJoints = self.holdSlideT.getJoints()
				#print "case_AB"
				self.holdT.reset(self.probeState, resultJoints)
		
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
					for i in range(self.numJoints-7, self.numJoints):				
						self.holdT.positions[i] = 180.0 / pi * stateJoints[i]
				else:
					for i in range(0, 6):					
						self.holdT.positions[i] = 180.0 / pi * stateJoints[i]
	
				self.holdT.step(self.probeState)
				
			else:
							
				if self.frontAnchoringState:
					
					if self.frontExtending:
						#self.doExtendFront()
						self.frontExtending = False
						self.frontExtendDone = True
		
						#if self.frontExtendDone:
						#	self.frontExtending = False
						
					else:
						#print "Front Anchor Step"

						self.doFrontAnchor()

						print "after doFrontAnchor, received frontAnchoringDone =", self.frontAnchoringDone
	
						if self.frontAnchoringDone:
							print "front anchoring done with splice joint =", self.spliceJoint
							self.frontAnchoringDone = False
							self.frontAnchoringState = False
							
				else:
					peakJoints = self.concertinaFit.getPeakJoints(self.currPeak)	
					
					print "Concertina Step"
					self.doBackConcertina()


				solids = []
				for i in range(self.numJoints):
					solids.append(self.concertinaFit.isJointSolid(i))
				#print "solids =", solids
				
				#print "joints =", self.spliceFitJoints()

				#print torques

				
				if not self.frontExtending:	
					self.transDone = self.holdT.step(self.probeState)
				else:
					self.transDone = self.holdSlideT.step(self.probeState)
		
				self.refDone = False
				
		else:
			
			" make no movements until all reference nodes are activated "
			if self.refDone:
				if not self.frontExtending:
					self.transDone = self.holdT.step(self.probeState)
				else:
					self.transDone = self.holdSlideT.step(self.probeState)

		" compute the mask "
		self.computeMaskAndOutput()
		
		if self.frontAnchoringState and self.frontCurve.isOutOfSegments(self.frontAnchorFit.lastPosition):

			" if we've exhausted too many joints, lets just quit with what we have "
			#if self.direction and self.spliceJoint >= 15:
			#	self.frontAnchoringDone = True
			#elif not self.direction and self.spliceJoint <= 23:
			#	self.frontAnchoringDone = True

			if self.lastAttempt:
				#self.frontAnchoringDone = True
				self.frontAnchoringDone = False
				self.frontAnchoringState = False

				print "ERROR: Quit with outOfSegments error on Front Anchoring"
				print "received frontAnchoringDone =", self.frontAnchoringDone
				
			else:
				
				if self.direction:
					self.spliceJoint += 2
					
					if self.spliceJoint >= self.topJoint:
						self.spliceJoint = self.topJoint
						self.lastAttempt = True

				else:
					self.spliceJoint -= 2
					if self.spliceJoint <= self.topJoint:
						self.spliceJoint = self.topJoint
						self.lastAttempt = True
					
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
				
				" compute the mask "
				#self.computeMaskAndOutput()
			
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
						
			self.lastAttempt = False
			
			if self.direction:
				self.spliceJoint = 7
			else:
				self.spliceJoint = 31
			
			"reset everything "
			self.frontAnchorFit = 0
			self.frontCurve = 0
			self.adaptiveCurve = 0
			self.concertinaFit = 0

			self.computeCurve()
			

		return isDone
	
	def computeMaskAndOutput(self):


		self.mask = [0.0 for i in range(0,self.numJoints)]

		anchorJoints = self.frontAnchorFit.getJoints()
		

		if self.frontAnchoringState and self.frontExtending:
			joints = self.holdSlideT.getJoints()
		else:
			joints = self.holdT.getJoints()
		
		if self.frontAnchoringState and self.frontExtendDone:
			" collect joint settings for both behaviors "
			#joints2 = self.globalCurveFit.getJoints()
			
			joints2 = [0.0 for i in range(self.numJoints)]
			
			if self.frontAnchorFit.anterior:			
				for p in range(0,self.frontAnchorFit.lastJoint+1):
					joints[p] = joints2[p]
			else:
				for p in range(self.frontAnchorFit.lastJoint, self.numJoints):
					joints[p] = joints2[p]

		" Set the joints "
		if self.refDone:

			self.mergeJoints([joints])
			#print "joints =", joints
			#print "frontAnchoringState:", self.frontAnchoringState
			
			" update mask for ContactReference "
			for i in range(self.numJoints):

				if (not self.frontAnchoringState and anchorJoints[i] != None):
					self.mask[i] = 1.0

				elif self.concertinaFit.isJointSolid(i):
					self.mask[i] = 1.0

				elif joints[i] == None:
					if self.frontAnchoringState:
						if self.direction and i > self.topJoint:
							self.mask[i] = 1.0
						
						if not self.direction and i < self.topJoint:
							self.mask[i] = 1.0

					else:
						self.mask[i] = 1.0			
					
		else:
			
			" transition the masks, and wait for them to be stable, then change the joints "
			
			allActive = True
			" update mask for ContactReference "
			for i in range(self.numJoints):
				if (not self.frontAnchoringState and anchorJoints[i] != None):
					self.mask[i] = 1.0
					
					" check if all reference points have been created "
					if not self.contacts.activeRef[i]:
						allActive = False
											
				elif self.concertinaFit.isJointSolid(i):
					self.mask[i] = 1.0

					" check if all reference points have been created "
					if not self.contacts.activeRef[i]:
						allActive = False

				elif joints[i] == None:
					if self.frontAnchoringState:
						if self.direction and i > self.topJoint:
							self.mask[i] = 1.0

							" check if all reference points have been created "
							if not self.contacts.activeRef[i]:
								allActive = False
						
						if not self.direction and i < self.topJoint:
							self.mask[i] = 1.0

							" check if all reference points have been created "
							if not self.contacts.activeRef[i]:
								allActive = False

					else:
						self.mask[i] = 1.0

						" check if all reference points have been created "
						if not self.contacts.activeRef[i]:
							allActive = False
							
			if allActive:
				self.refDone = True


	def reset(self, probeState):

		self.probeState = probeState

		self.mask = [0.0 for i in range(0,39
									)]
		self.frontAnchoringState = True
		self.frontAnchoringDone = False
		self.frontExtending = True
		self.frontExtendDone = False
					
		self.minAmp = 0.0
		self.maxAmp = 0.0			
		self.ampInc = 0.04
					
		self.currPeak = 0
		self.lastSpliceAngle = 0.0
					
		self.lastAttempt = False
		
		if self.direction:
			self.spliceJoint = 7
		else:
			self.spliceJoint = 31
		
		"reset everything "
		self.frontAnchorFit = 0
		self.frontCurve = 0
		self.adaptiveCurve = 0
		self.concertinaFit = 0

		self.computeCurve()
			
		#self.computeCenterPoints()
		#self.frontCurve.saveAmps()
		
		"reset everything "
		
		#self.frontAnchorFit = 0
		#self.frontCurve = 0
		#self.computeCurve()
		
		#self.mask = [0.0 for i in range(0,40)]
		#self.count = 0
	
		self.currPeak = 0
	
	
	
	
