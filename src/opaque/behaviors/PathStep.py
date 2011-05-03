

from Behavior import Behavior
from VoronoiFit import VoronoiFit
from FrontAnchorFit import FrontAnchorFit
from FastLocalCurveFit import FastLocalCurveFit
from HoldTransition import HoldTransition
from HoldSlideTransition import HoldSlideTransition
from AdaptiveAnchorCurve import AdaptiveAnchorCurve
from BackConcertinaCurve import BackConcertinaCurve
from GlobalCurveFit import GlobalCurveFit
from GlobalCurveSlide import GlobalCurveSlide

from copy import copy, deepcopy
from math import pi, cos, sin, fabs, sqrt

"""
1.  Check condition in FrontAnchorFit when curve goes beyond joint 0
2.  If joint 0 passed and anchor is loose, start over with higher splice joint
3.  Changing splice joint will shorten tail of BackConcertinaCurve
4.  Add the front anchoring state and then the concertina gait state to behavior control
"""

class PathStep(Behavior):

	def __init__(self, robotParam, probeState, contacts, mapGraph, direction = True, path = []):
		Behavior.__init__(self, robotParam)

		print "creating PathStep"

		self.numJoints = self.robotParam['numJoints']
		self.numSegs = self.robotParam['numSegs']
		self.maxTorque = self.robotParam['maxTorque']
		self.probeState = probeState

		self.mapGraph = mapGraph
		#self.localNode = self.mapGraph.getCurrentNode()
		self.contacts = contacts
		self.frontAnchorFit = 0
		self.concertinaFit = 0

		self.compliantTorque = 0.005
		#self.compliantTorque = 3.0
		
		self.direction = direction
		self.isInit = False

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
		self.jerkState = 0

		self.errorJoint = 0
		self.errorPose1 = [0.0,0.0,0.0]
		self.errorPose2 = [0.0,0.0,0.0]
		self.originPose = [0.0,0.0,0.0]
		
		#self.jerkState = -1
		self.nomJerk = 0.0
		self.nomJerks = [0.0,0.0,0.0,0.0]
		self.nomJerks2 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
		self.jerkAngles = [0.0,0.0,0.0,0.0]
		self.jerkAngles2 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
		self.prevJerkAngles = [0.0,0.0,0.0,0.0]
		
		
		#self.nomJerk = 0.0
		#self.nomJerks = [0.0,0.0,0.0,0.0]
		#self.jerkAngles = [0.0,0.0,0.0,0.0]
		#self.prevJerkAngles = [0.0,0.0,0.0,0.0]
		
		if direction:
			self.spliceJoint = 7
		else:
			self.spliceJoint = 31
			
		self.lastSpliceAngle = 0.0

		self.mask = [0.0 for i in range(0,40)]
		self.count = 0
		
		self.frozenJoints = [None for i in range(self.numJoints)]
		
		self.currPeak = 0
		
		self.minAmp = 0.0
		self.maxAmp = 0.0
		self.ampInc = 0.04
			
		self.plotCount = 0
		
		" flag true when the HoldTransition behavior has completed its behavior, select next transition "
		self.transDone = False

		" flag true when all the joints are stable, then step the HoldTransition behavior "
		self.refDone = False

		self.frontAnchoringState = True
		self.frontAnchoringDone = False
		self.frontExtending = True
		self.frontExtendDone = False
				
	def getMask(self):
		return self.mask

	def setPath(self, path):
		print "setting path: ", path

		self.path = deepcopy(path)
		#self.computeCurve()

		self.pathCurve1 = VoronoiFit(deepcopy(self.path))
		self.pathCurve2 = VoronoiFit(deepcopy(self.path))
		
		print "creating GlobalCurveFit"
		#self.globalCurveFit = GlobalCurveFit(self.robotParam, self.contacts, self.pathCurve1, localNode = self.mapGraph.currNode)
		self.globalCurveFit = GlobalCurveFit(self.robotParam, self.contacts, self.pathCurve1)
		#self.globalCurveFit.setTimerAliasing(10)
		
		print "creating GlobalCurveSlide"
		#self.globalCurveSlide = GlobalCurveSlide(self.robotParam, self.contacts, self.pathCurve2, localNode = self.mapGraph.currNode)
		self.globalCurveSlide = GlobalCurveSlide(self.robotParam, self.contacts, self.pathCurve2)
		
		#self.setDirection(not self.globalCurveFit.getPathDirection())
		result = self.globalCurveSlide.getPathDirection()
		
		print "received", result, "from globalCurveSlide"
		
		self.setDirection(result)
		
		#self.globalCurveFit.clearDraw()
		#self.globalCurveSlide.clearDraw()
		
	def reverseDirection(self):
		self.path.reverse()
		
		self.setPath(self.path)
		
		self.computeCurve()
						
	def setDirection(self, isForward):

		self.direction = isForward
		
		print "setting PathStep direction =", self.direction

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
				
	def computeCurve(self):

				
		self.frontCurve = AdaptiveAnchorCurve(4*pi, self.robotParam['segLength'])				

		if self.frontAnchorFit == 0:
			self.frontAnchorFit = FrontAnchorFit(self.robotParam, self.direction, self.spliceJoint)
			self.frontAnchorFit.setCurve(self.frontCurve)
			self.frontAnchorFit.setSpliceJoint(self.spliceJoint)

		else:
			self.frontAnchorFit.setCurve(self.frontCurve)
			self.frontAnchorFit.setSpliceJoint(self.spliceJoint)


		self.adaptiveCurve = BackConcertinaCurve(4*pi)
		self.adaptiveCurve.setTailLength(2.0)
		
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

		self.frontAnchorFit.step(self.probeState)
		self.concertinaFit.step(self.probeState)
		
		resultJoints = self.spliceFitJoints()
		
		" straighten the part of frontAnchorFit that is not actuated "
		if self.frontAnchorFit.anterior:
			startNode = 0
			endNode = self.frontAnchorFit.lastJoint			
		else:
			startNode = self.frontAnchorFit.lastJoint
			endNode = self.numJoints-1
						
		self.globalCurveFit.setBoundaries(startNode, endNode)
			
		for i in range(startNode, endNode+1):
			#print "changing joint", i, "from", resultJoints[i], "to 0.0"
			resultJoints[i] = 0.0
			

		
		#self.holdSlideT.setSpliceJoint(self.spliceJoint)

		self.globalCurveSlide.reset(self.probeState, resultJoints)
		
		self.globalCurveSlide.step(self.probeState)
		resultJoints2 = self.globalCurveSlide.getJoints()
		
		self.globalCurveFit.step(self.probeState)
		
		self.holdT.reset(self.probeState, resultJoints2)
		self.refDone = False
		
		#self.globalCurveFit.clearDraw()
		#self.globalCurveSlide.clearDraw()
		
		self.computeMaskAndOutput()

	def spliceFitJoints(self, isStart = False):
		
		result = [None for i in range(self.numJoints)]
		
		if self.frontAnchoringState:
			joints1 = self.frontAnchorFit.getJoints()
		else:
			joints1 = self.frozenJoints
		
		joints2 = self.concertinaFit.getJoints()
		
		self.concertinaFit.spliceJoint = self.spliceJoint
		self.concertinaFit.spliceAngle = pi / 180.0 * joints1[self.spliceJoint]
				
		newJointVal = joints1[self.spliceJoint] + joints2[self.spliceJoint]
				
		tempJoints = copy(joints1)
		tempJoints[self.spliceJoint] = newJointVal
		
		#print "splicing:"
		#print tempJoints
		#print joints2
		
		jointList = [tempJoints,joints2]
		
		" merge the joints prioritizing the first to last "
		for joints in jointList:
			for i in range(len(joints)):
				if result[i] == None and joints[i] != None:
					result[i] = joints[i]
		
		return result		

	"""
	def computeCmdSegPoints(self):

		segLength = self.robotParam['segLength']

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
				xTotal = cmdOrigin[0] - segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] - segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.insert(0, cmdOrigin)

			cmdOrigin = [segLength*cos(-spliceJointAngle),segLength*sin(-spliceJointAngle),-spliceJointAngle]
			self.cmdSegPoints.append(cmdOrigin)

			for i in range(self.spliceJoint+1,self.numJoints):
				
				sampAngle = self.probe.getServoCmd(i)

				totalAngle = cmdOrigin[2] - sampAngle
				xTotal = cmdOrigin[0] + segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] + segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.append(cmdOrigin)
			
		else:
			
			spliceJointAngle = -self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0 
			
			self.cmdSegPoints = []
			cmdOrigin = [segLength*cos(-spliceJointAngle),segLength*sin(-spliceJointAngle),-spliceJointAngle + pi]

			self.cmdSegPoints.append(cmdOrigin)

			ind = range(0,self.spliceJoint)
			ind.reverse()
			
			for i in ind:
				sampAngle = self.probe.getServoCmd(i)
				totalAngle = cmdOrigin[2] + sampAngle
				xTotal = cmdOrigin[0] - segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] - segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.insert(0,cmdOrigin)

			cmdOrigin = [0.0,0.0,-spliceJointAngle + pi]			
			self.cmdSegPoints.append(cmdOrigin)
			
			for i in range(self.spliceJoint,self.numJoints):
				sampAngle = self.probe.getServoCmd(i)
				totalAngle = cmdOrigin[2] - sampAngle
				xTotal = cmdOrigin[0] + segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] + segLength*sin(totalAngle)
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

		segLength = self.robotParam['segLength']

		if self.direction:
			" the actual configuration of the snake "
			#actSegPoints = [[-segLength,0.0,0.0]]

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
				xTotal = actOrigin[0] - segLength*cos(totalAngle)
				zTotal = actOrigin[1] - segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.insert(0,actOrigin)
	
			actOrigin = [segLength*cos(-spliceJointAngle),segLength*sin(-spliceJointAngle),-spliceJointAngle]
			actSegPoints.append(actOrigin)
			
			for i in range(self.spliceJoint+1,self.numJoints):
				sampAngle = stateJoints[i]
				totalAngle = actOrigin[2] - sampAngle
				xTotal = actOrigin[0] + segLength*cos(totalAngle)
				zTotal = actOrigin[1] + segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.append(actOrigin)

		else:
			
			spliceJointAngle = -self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0 
			
			actSegPoints = []
			actOrigin = [segLength*cos(-spliceJointAngle),segLength*sin(-spliceJointAngle),-spliceJointAngle + pi]

			actSegPoints.insert(0, actOrigin)

			ind = range(0,self.spliceJoint)
			ind.reverse()
			
			for i in ind:
				sampAngle = stateJoints[i]
				totalAngle = actOrigin[2] + sampAngle
				xTotal = actOrigin[0] - segLength*cos(totalAngle)
				zTotal = actOrigin[1] - segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.insert(0,actOrigin)

			actOrigin = [0.0,0.0,-spliceJointAngle + pi]			
			actSegPoints.append(actOrigin)
			
			for i in range(self.spliceJoint,self.numJoints):
				sampAngle = stateJoints[i]
				totalAngle = actOrigin[2] - sampAngle
				xTotal = actOrigin[0] + segLength*cos(totalAngle)
				zTotal = actOrigin[1] + segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.append(actOrigin)
				
		if self.direction:

			
			spliceJointAngle = self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0

			" local peak configurations "
			peakCurves = [[]]

			peakOrigin = [segLength*cos(-spliceJointAngle),segLength*sin(-spliceJointAngle),-spliceJointAngle]
			peakCurves[-1].append(copy(peakOrigin))


			for i in range(self.spliceJoint+1,self.numJoints):
				sampAngle = stateJoints[i]
				totalAngle = peakOrigin[2] - sampAngle
				xTotal = peakOrigin[0] + segLength*cos(totalAngle)
				zTotal = peakOrigin[1] + segLength*sin(totalAngle)
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
					xTotal = peakOrigin[0] + segLength*cos(totalAngle)
					zTotal = peakOrigin[1] + segLength*sin(totalAngle)
					pnt = [xTotal, zTotal, totalAngle]
					peakOrigin = pnt
					peakCurves[-1].append(copy(peakOrigin))
				else:
					peakCurves[-1].append(copy(peakOrigin))
		else:

			spliceJointAngle = -self.concertinaFit.getJoints()[self.spliceJoint] * pi /180.0

			" local peak configurations "
			peakCurves = [[]]

			peakOrigin = [segLength*cos(-spliceJointAngle),segLength*sin(-spliceJointAngle),-spliceJointAngle + pi]
			peakCurves[-1].append(copy(peakOrigin))
			
			ind = range(0,self.spliceJoint)
			ind.reverse()

			for i in ind:
				sampAngle = stateJoints[i]
				totalAngle = peakOrigin[2] + sampAngle
				xTotal = peakOrigin[0] - segLength*cos(totalAngle)
				zTotal = peakOrigin[1] - segLength*sin(totalAngle)
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
					xTotal = peakOrigin[0] - segLength*cos(totalAngle)
					zTotal = peakOrigin[1] - segLength*sin(totalAngle)
					pnt = [xTotal, zTotal, totalAngle]
					peakOrigin = pnt
					peakCurves[-1].append(copy(peakOrigin))
				else:
					peakCurves[-1].append(copy(peakOrigin))
					
			
		errors = []
		for i in range(0,self.probe.numSegs-2):
			errors.append(fabs(stateJoints[i]-self.probe.getServoCmd(i)))

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
		

		#print "plotting graph", self.plotCount

		pylab.xlim(-4,4)
		#pylab.xlim(-0.5,4)
		pylab.ylim(-3,3)
		pylab.savefig("fitPlot%04u.png" % self.plotCount)
		
		self.plotCount += 1
	"""	
		
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

		err1 = sqrt(self.errorPose1[0]**2 + self.errorPose1[1]**2)

		numPoses = len(self.errorPoses2)
		xDiff = 0.0
		yDiff = 0.0
		for p in self.errorPoses2:
			xDiff += p[0]
			yDiff += p[1]
			
		xDiff /= numPoses
		yDiff /= numPoses
		
		err2 = sqrt(xDiff**2 + yDiff**2)

		err2 = sqrt(self.errorPose2[0]**2 + self.errorPose2[1]**2)

		numPoses = len(self.errorPoses3)
		xDiff = 0.0
		yDiff = 0.0
		for p in self.errorPoses3:
			xDiff += p[0]
			yDiff += p[1]
			
		xDiff /= numPoses
		yDiff /= numPoses
		
		err3 = sqrt(xDiff**2 + yDiff**2)
				
		return err1, err2, err3
			
	def doExtendFront(self):

		print "Extend Front"

		if self.globalCurveSlide.step(self.probeState):
			self.frontExtendDone = True
			self.frontExtending = False
			#self.globalCurveSlide.clearDraw()
				
		joints = self.globalCurveSlide.getJoints()

		self.holdT.reset(self.probeState, joints)
						
	def doFrontAnchor(self):
		
		stateJoints = self.probeState['joints']
		cmdJoints = self.probeState['cmdJoints']
		errors = self.probeState['errors']

		" compute the maximum error of the joints in the current peak "
		anchorJoints = self.frontAnchorFit.getPeakJoints()
		anchorJoints.sort()

		maxError = 0.0	
		for j in anchorJoints:
			err = errors[j]
			if fabs(err) > maxError:
				maxError = fabs(err)
		
		" check for terminating criteria "
		currAmp = self.frontCurve.getPeakAmp()	
		nextVal = currAmp	

		" draw the state of the curve fitting "
		#self.drawFit()
		
		if not self.isJerking:
			" perform amplitude operation here "

			#print "amps before:", self.minAmp, self.maxAmp, self.ampInc, currAmp, nextVal
			
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
			
			#print "amps after:", self.minAmp, self.maxAmp, self.ampInc, currAmp, nextVal
				
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
						
					self.errorJoint = maxJoint

					"""
					self.jerkJoints2 = []
					for k in range(4,12):
						jointNum = maxJoint + k + 1
						if jointNum > 38:
							break
						self.jerkJoints2.append(jointNum)
					"""

				else:
					minJoint = min(anchorJoints)
					self.jerkJoint = minJoint - 1
					
					self.jerkJoints = []
					for k in range(0,4):
						jointNum = minJoint - k - 1
						if jointNum < 0:
							break
						self.jerkJoints.append(jointNum)

					self.errorJoint = minJoint

					"""
					self.jerkJoints2 = []
					for k in range(4,12):
						jointNum = minJoint - k - 1
						if jointNum < 0:
							break
						self.jerkJoints2.append(jointNum)
					"""

				#print "anchor joints =", anchorJoints
				#print "setting jerkJoints to", self.jerkJoints
				
				if self.jerkState == -1:
					
					for k in range(len(self.jerkJoints2)):
						self.nomJerks2[k] = stateJoints[self.jerkJoints2[k]] * 180.0 / pi
				
					for k in range(len(self.jerkJoints2)):
						if k % 2 == 0:
							self.jerkAngles2[k] = 80
						else:
							self.jerkAngles2[k] = -80

					self.jerkState = 0
				
				elif self.jerkState == 0:
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

					self.originPose = self.contacts.getAveragePose(self.errorJoint)

					#self.originGndPoses = []
					#for k in anchorJoints:
					#	self.originGndPoses.append(self.contacs.getActualJointPose(k))

					self.jerkState = 1

				elif self.jerkState == 1:
					
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

					tempPose = self.contacts.getAveragePose(self.errorJoint)
					self.errorPose1 = [tempPose[0]-self.originPose[0],tempPose[1]-self.originPose[1],tempPose[2]-self.originPose[2]]
				
					#print self.errorPoses1
					
					self.originPoses = []
					for k in anchorJoints:
						self.originPoses.append(self.contacts.getAveragePose(k))

					self.jerkState = 2
					
				elif self.jerkState == 2:
					#print "jerk angle error = " , stateJoints[self.jerkJoint]-cmdJoints[self.jerkJoint]
					self.jerkErrors.append(stateJoints[self.jerkJoint]-cmdJoints[self.jerkJoint])

					" error = (-0.09, 0.005) is good "
					" error = (-2.73, -1.99) is bad "
					self.prevJerkAngle = self.jerkAngle
					self.jerkAngle = 0
					
					for k in range(len(self.jerkJoints)):
						self.prevJerkAngles[k] = self.jerkAngles[k]
						self.jerkAngles[k] = 0

					#for k in range(len(self.jerkJoints2)):
					#	self.jerkAngles2[k] = 0

					self.errorPoses2 = []
					for k in range(len(anchorJoints)):
						tempPose = self.contacts.getAveragePose(anchorJoints[k])
						originPose = self.originPoses[k]
						self.errorPoses2.append([tempPose[0]-originPose[0],tempPose[1]-originPose[1],tempPose[2]-originPose[2]])

					tempPose = self.contacts.getAveragePose(self.errorJoint)
					self.errorPose2 = [tempPose[0]-self.originPose[0],tempPose[1]-self.originPose[1],tempPose[2]-self.originPose[2]]

					#self.jerkState = 3
					self.jerkState = 5
					self.errorPoses3 = [[0.0,0.0,0.0] for i in range(4)]
				
					#print self.errorPoses2

				elif self.jerkState == 3:
					self.prevJerkAngle = self.jerkAngle

					for k in range(len(self.jerkJoints)):
						self.prevJerkAngles[k] = self.jerkAngles[k]

					for k in range(len(self.jerkJoints)):
						if k % 2 == 0:
							self.jerkAngles[k] = 80
						else:
							self.jerkAngles[k] = -80
						
					self.jerkState = 4
				
				elif self.jerkState == 4:
					self.prevJerkAngle = self.jerkAngle

					for k in range(len(self.jerkJoints)):
						self.prevJerkAngles[k] = self.jerkAngles[k]
						self.jerkAngles[k] = 0

					self.errorPoses3 = []
					for k in range(len(anchorJoints)):
						tempPose = self.contacts.getAveragePose(anchorJoints[k])
						originPose = self.originPoses[k]
						self.errorPoses3.append([tempPose[0]-originPose[0],tempPose[1]-originPose[1],tempPose[2]-originPose[2]])

					#print "errorPoses3:", self.errorPoses3

					self.jerkState = 5
										
				elif self.jerkState == 5:
					self.jerkingDone = True

					for k in range(len(self.jerkJoints)):
						self.prevJerkAngles[k] = self.jerkAngles[k]
						
					self.jerkState = 0
					#self.jerkState = -1

				#print "setting jerk joint to ", self.jerkAngle + self.nomJerk
				#print "jerk angle error = " , self.probe.getServo(j)-self.probe.getServoCmd(j)

		if self.isJerking and self.jerkingDone:
			
			" if the anchor is not secure, lets try it again "
			err1, err2, err3 = self.computeAnchorErrors()
			
			#print "errorPoses1:", self.errorPoses1
			#print "errorPoses2:", self.errorPoses2
			#print "errorPoses3:", self.errorPoses3
			print "anchor errors =", err1, err2, err3
			
			#if err1 > 0.1 or err2 > 0.1:
			if False:

				" TODO " 
				" if the anchor failed, increment the splice joint higher "

				#if self.direction:
				#	self.spliceJoint = 7
				#else:
				#	self.spliceJoint = 31
				
				" if we've exhausted too many joints, lets just quit with what we have "
				if self.direction and self.spliceJoint >= 15:
					self.frontAnchoringDone = True
				elif not self.direction and self.spliceJoint <= 23:
					self.frontAnchoringDone = True
				
				else:
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

		" Lead joints of the snake are weakened when not part of front-anchoring "		
		if len(anchorJoints) > 0 and self.frontCurve.getPeakAmp() > 0:
			
			if not self.direction:
				maxJoint = max(anchorJoints)
				#print "weakening joints", range(maxJoint+1, self.numJoints)
				for i in range(maxJoint+1, self.numJoints):
					if i <= self.numJoints-1:
						self.torques[i] = self.compliantTorque
			else:
				minJoint = min(anchorJoints)
				#print "weakening joints", range(0,minJoint-1)
				for i in range(0, minJoint-1):					
					if i <= self.numJoints-1:
						self.torques[i] = self.compliantTorque
	
		" execute the local curve fitting "
		self.frontAnchorFit.step(self.probeState)

		" compute the bounds for behavior for global curve fitting "
		if self.frontAnchorFit.anterior:
			startNode = 0
			endNode = self.frontAnchorFit.lastJoint
			self.globalCurveFit.setBoundaries(startNode, endNode)
		else:
			startNode = self.frontAnchorFit.lastJoint
			endNode = self.numJoints-1	
			self.globalCurveFit.setBoundaries(startNode, endNode)

		
		#self.globalCurveFit.step()
			
		" Set the joints, with joints2 as priority over joints1 "
		#self.mergeJoints([joints2, joints1])

		resultJoints = self.spliceFitJoints()
	
		self.holdT.reset(self.probeState, resultJoints)
		
		if self.isJerking:
			for k in range(len(self.jerkJoints)):
				self.holdT.positions[self.jerkJoints[k]] = self.nomJerks[k] + self.jerkAngles[k]
			
	def doBackConcertina(self):
		
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
		
		currAmp = self.adaptiveCurve.getPeakAmp(self.currPeak)
		nextVal = currAmp

		" perform amplitude operation here "
		#print "errors:", maxError, errors
		#print self.probeState['joints']
		#print self.probeState['cmdJoints']
		#print "amps before:", self.minAmp, self.maxAmp, self.ampInc, currAmp, nextVal
		
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

		#print "amps after:", self.minAmp, self.maxAmp, self.ampInc, currAmp, nextVal

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
						self.torques[i] = self.compliantTorque
			else:
				minJoint = min(peakJoints)

				for i in range(0, minJoint-1):					
					if i <= self.numJoints-1:
						self.torques[i] = self.compliantTorque

		resultJoints = self.spliceFitJoints()
		self.holdT.reset(self.probeState, resultJoints)

		peakJoints = self.concertinaFit.getPeakJoints(self.currPeak) + self.concertinaFit.getPeakJoints(self.currPeak+1)
		peakJoints.sort()
		
	def step(self, probeState):
		Behavior.step(self, probeState)
				
		self.count += 1
		isDone = False

		stateJoints = self.probeState['joints']		
		torques = self.probeState['torques']

		self.torques = []
		for i in range(0,self.numJoints):
			self.torques.append(torques[i])
		
		if self.transDone:
			self.transDone = False

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
					
				joints = probeState['joints']
				
				" prevent the back anchor from slithering around or buckling "
				if self.direction:
					for i in range(self.numJoints-9, self.numJoints-1):				
						self.holdT.positions[i] = 180.0 / pi * joints[i]
				else:
					for i in range(0, 6):					
						self.holdT.positions[i] = 180.0 / pi * joints[i]
	
				self.holdT.step(probeState)
				
			else:
							
				if self.frontAnchoringState:
					
					if self.frontExtending:
						
						print "Front Extend Step"
						self.doExtendFront()
						if self.frontExtendDone:
							self.frontExtending = False
						
					else:
						print "Front Anchor Step"
						self.doFrontAnchor()
		
						if self.frontAnchoringDone:
							self.frontAnchoringDone = False
							self.frontAnchoringState = False
							self.frontExtended = False
							
							" save the current joints and use them from now on "

							joints = probeState['joints']
							
							" prevent the front anchor from causing error "
							if self.direction:
								for i in range(0, self.spliceJoint):				
									self.frozenJoints[i] = 180.0 / pi * joints[i]
								self.frozenJoints[self.spliceJoint] = self.frontAnchorFit.getJoints()[self.spliceJoint]
							else:
								for i in range(self.spliceJoint+1, self.numJoints):					
									self.frozenJoints[i] = 180.0 / pi * joints[i]
								self.frozenJoints[self.spliceJoint] = self.frontAnchorFit.getJoints()[self.spliceJoint]
					
				else:
					peakJoints = self.concertinaFit.getPeakJoints(self.currPeak)	
					
					print "Concertina Step"
					self.doBackConcertina()


				self.refDone = False
				self.transDone = self.holdT.step(probeState)
				
		else:
			
			" make no movements until all reference nodes are activated "
			if self.refDone:
				self.transDone = self.holdT.step(probeState)					
				
		self.computeMaskAndOutput()

		
		if self.frontCurve.isOutOfSegments(self.frontAnchorFit.lastPosition):

			print "NOTE: ran out of segments on splice joint", self.spliceJoint
			
			" if we've exhausted too many joints, lets just quit with what we have "
			if self.direction and self.spliceJoint >= 15:
				self.frontAnchoringDone = True
			elif not self.direction and self.spliceJoint <= 23:
				self.frontAnchoringDone = True
			else:
				
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
			

		return isDone
	
	def computeMaskAndOutput(self):
		
		anchorJoints = self.frontAnchorFit.getJoints()
		
		joints = self.holdT.getJoints()
		
		if self.frontAnchoringState and self.frontExtendDone:
			
			self.globalCurveFit.step(self.probeState)
			
			" collect joint settings for both behaviors "
			
			joints2 = self.globalCurveFit.getJoints()			
			
			if self.frontAnchorFit.anterior:			
				for p in range(0,self.frontAnchorFit.lastJoint+1):
					joints[p] = joints2[p]
			else:
				for p in range(self.frontAnchorFit.lastJoint, self.numJoints):
					joints[p] = joints2[p]

		segLength = self.robotParam['segLength']

		" Set the joints "
		if self.refDone:

			self.mergeJoints([joints])
			" update mask for ContactReference "
			for i in range(self.numJoints):

				if (not self.frontAnchoringState and anchorJoints[i] != None) or self.concertinaFit.isJointSolid(i) or joints[i] == None:
					self.mask[i] = 1.0
				
					" for activating the back anchor joints when front anchoring or extending "
				#elif self.frontAnchoringState and self.direction and i > (self.spliceJoint + 2.0/segLength):
				elif self.frontAnchoringState and self.direction and i > 19:
					
					self.mask[i] = 1.0
							
				#elif self.frontAnchoringState and not self.direction and i < (self.spliceJoint - 2.0/segLength):
				elif self.frontAnchoringState and not self.direction and i < 19:
					self.mask[i] = 1.0
				
				else:
					self.mask[i] = 0.0
								
		else:
			" transition the masks, and wait for them to be stable, then change the joints "

			allActive = True
			" update mask for ContactReference "
			for i in range(self.numJoints):
				if (not self.frontAnchoringState and anchorJoints[i] != None) or self.concertinaFit.isJointSolid(i) or joints[i] == None:
					self.mask[i] = 1.0

					" check if all reference points have been created "
					if not self.contacts.activeRef[i]:
						allActive = False
				
				" for activating the back anchor joints when front anchoring or extending "
				if self.frontAnchoringState:
					#if self.direction and i > (self.spliceJoint + 2.0/segLength):
					if self.direction and i > 19:
						self.mask[i] = 1.0

						" check if all reference points have been created "
						if not self.contacts.activeRef[i]:
							allActive = False
							
						" add the mask for the inactive portion of the back concertina gait while front anchoring "
					#elif not self.direction and i < (self.spliceJoint - 2.0/segLength):
					elif not self.direction and i < 19:
						self.mask[i] = 1.0
				
						" check if all reference points have been created "
						if not self.contacts.activeRef[i]:
							allActive = False
					
			if allActive:
				self.refDone = True

	def reset(self, probeState):
		
		self.probeState = probeState

		"reset everything "
		
		self.frontAnchorFit = 0
		self.frontCurve = 0
		self.computeCurve()
		
		self.mask = [0.0 for i in range(0,40)]
		self.count = 0
	
		self.currPeak = 0
	
	
	
	
