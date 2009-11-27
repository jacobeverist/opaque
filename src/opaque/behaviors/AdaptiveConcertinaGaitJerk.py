import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from copy import *
from common import *
from Behavior import *
from FastLocalCurveFit import FastLocalCurveFit
from HoldTransition import HoldTransition

class AdaptiveConcertinaGaitJerk(Behavior):

	def __init__(self, probe, contacts, mapGraph, direction = True):
		Behavior.__init__(self, probe)

		self.mapGraph = mapGraph
		self.localNode = self.mapGraph.getCurrentNode()
		self.contacts = contacts
		self.localFit = 0
		self.direction = direction
		self.isInit = False

		self.holdT = HoldTransition(probe)

		self.cPoints = []

		self.minWidth = 0.0
		self.maxWidth = 0.0
		self.widthThresh = 0.0
		self.incWidth = False
		self.decWidth = False
		self.nomWidth = 0.0
		self.widthInc = 0.04
		self.isWidthMode = False
		self.widthLeft = False
		
		self.isJerking = False
		self.jerkingDone = False
		self.jerkAngle = 0
		self.prevJerkAngle = self.jerkAngle
		self.jerkErrors = []
		
		self.computeCurve()

		self.mask = [0.0 for i in range(0,40)]
		self.count = 0
		
		self.currPeak = 0
		
		self.minAmp = 0.0
		self.maxAmp = 0.0

		self.ampInc = 0.04
		

		
		self.isDone = False
		self.plotCount = 0
		
		self.transDone = False
		self.newJoints = [None for i in range(0,self.probe.numSegs-1)]
		self.refDone = False
		
	def getMask(self):
		return self.mask
		
	def setDirection(self, isForward):

		if isForward:
			self.localFit.setSide(True)
			self.localFit.setRootNode(38)
		else:
			self.localFit.setSide(False)
			self.localFit.setRootNode(0)
		
	def computeCurve(self):
		self.curve = AdaptiveCosine(4*pi, 1.0)
		#self.curve.setHeadLength(0.5)
		#self.curve = AdaptiveCosine(2*pi, 1.0)
		
		if self.localFit == 0:
			self.localFit = FastLocalCurveFit(self.probe, self.direction, self.curve, 38)
		else:
			self.localFit.setCurve(self.curve)

		self.holdT.reset(self.localFit.getJoints())
		self.refDone = False

		self.nomWidth = self.curve.ampSpacing
		self.minWidth = self.nomWidth
		self.maxWidth = self.nomWidth
		self.widthThresh = self.nomWidth + 3*self.widthInc
			
	def getCenterPoints(self):
		
		#return copy(self.cmdSegPoints)
		return self.cPoints
		
	
	def computeCenterPoints(self):
		
		points = [[self.curve.ampSpacing/2.0 + i * self.curve.ampSpacing, -self.curve.peakAmp[0]] for i in range(self.curve.numPeaks)]

		" 1. at each inflection point, find the closest segment and it's point on segment"
		" 2. cycle through all 40 segments "
		" 3. transform those points to current estimated positions "

		self.localNode = self.mapGraph.getCurrentNode()
		
		c_configs = []
		
		#print len(points), "curve points:"
		#print points
		
		for p in points:
			minDist = 1e100
			minPoint = [0.0,0.0]
			jointIndex = 0
			segLen = 0.0
			
			#for i in range(len(self.solidJointPositions)-1):
			for i in range(len(self.cmdSegPoints)-1):
				#for j in [[-ARMLENGTH,0.0]] + self.solidJointPositions:
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

		#print len(c_configs), "configurations:"
		#print c_configs
		
		" return joint number and position on segment "
		#c_configs =  self.localFit.getCenterPoints()
		
		self.cPoints = []
		for c in c_configs:
			joint = c[0]
			segLen = c[1]
			
			#pose = self.contacts.getClosestPose(joint)
			pose = self.localNode.getJointPose(joint)
			point = [pose[0] + segLen*cos(pose[2]), pose[1] + segLen*sin(pose[2])]
			
			if point not in self.cPoints:
				self.cPoints.append(point)
		
		#print "returning", len(self.cPoints), "center points: "
		#print self.cPoints
	
	def computeCmdSegPoints(self):	
		
		if self.direction:
				
			" the configuration of the snake as we've commanded it "
			self.cmdSegPoints = [[-self.probe.segLength,0.0,0.0]]
			cmdOrigin = [0.0,0.0,0.0]
			self.cmdSegPoints.append(cmdOrigin)
	
			for i in range(0,self.probe.numSegs-1):
				
				if self.localFit.isJointSolid(i):
					sampAngle = self.localFit.getSolidJoint(i)
				else:
					sampAngle = self.probe.getServoCmd(i)

				totalAngle = cmdOrigin[2] - sampAngle
				xTotal = cmdOrigin[0] + self.probe.segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] + self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.append(cmdOrigin)

		else:
				
			" the configuration of the snake as we've commanded it "
			self.cmdSegPoints = [[-self.probe.segLength,0.0,0.0]]
			cmdOrigin = [0.0,0.0,pi]
			self.cmdSegPoints.insert(0, cmdOrigin)
		
			ind = range(0,self.probe.numSegs-1)
			ind.reverse()
	
			for i in ind:
				if self.localFit.isJointSolid(i):
					sampAngle = self.localFit.getSolidJoint(i)
				else:
					sampAngle = self.probe.getServoCmd(i)
				#sampAngle = self.probe.getServoCmd(i)

				totalAngle = cmdOrigin[2] + sampAngle
				xTotal = cmdOrigin[0] - self.probe.segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] - self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.insert(0, cmdOrigin)
				#cmdSegPoints[i+2] = cmdOrigin
			
	def drawFit(self):
		
		#peakJoints = self.localFit.getPeakJoints(self.currPeak)

		pylab.clf()
		crvPoints = self.curve.getPoints()		
		xP = []
		yP = []
		for p in crvPoints:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP,color='0.3')
		
		self.computeCmdSegPoints()

		if self.direction:
			" the actual configuration of the snake "
			actSegPoints = [[-self.probe.segLength,0.0,0.0]]
			actOrigin = [0.0,0.0,0.0]
			actSegPoints.append(actOrigin)
			
			for i in range(0,self.probe.numSegs-1):
				sampAngle = self.probe.getServo(i)
				totalAngle = actOrigin[2] - sampAngle
				xTotal = actOrigin[0] + self.probe.segLength*cos(totalAngle)
				zTotal = actOrigin[1] + self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.append(actOrigin)	
		else:
			" the actual configuration of the snake "
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
		

		if self.direction:
			" local peak configurations "
			peakCurves = [[]]
			peakCurves[-1].append([-self.probe.segLength,0.0,0.0])
			peakOrigin = [0.0,0.0,0.0]
			peakCurves[-1].append(copy(peakOrigin))

			for i in range(0,self.probe.numSegs-1):
				sampAngle = self.probe.getServo(i)
				totalAngle = peakOrigin[2] - sampAngle
				xTotal = peakOrigin[0] + self.probe.segLength*cos(totalAngle)
				zTotal = peakOrigin[1] + self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				peakOrigin = pnt
				#peakCurves[-1].append(copy(peakOrigin))
				
				" if current joint is the minimum of a current peak, reset the origin "
				peakIndex = -1
				for index in range(len(self.localFit.jointClasses)):
					if self.localFit.jointClasses[index].count(i) > 0:
						peakIndex = index
						break
				
				if peakIndex >= 1 and min(self.localFit.jointClasses[peakIndex]) == i:
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
						
					#peakOrigin = cmdSegPoints[2+i]
					#peakCurves[-1].append(copy(peakOrigin))
		else:
			" local peak configurations "
			peakCurves = [[]]
			peakCurves[-1].append([-self.probe.segLength,0.0,0.0])
			peakOrigin = [0.0,0.0,pi]
			peakCurves[-1].append(copy(peakOrigin))

			ind = range(0,self.probe.numSegs-1)
			ind.reverse()
			
			for i in ind:
				sampAngle = self.probe.getServo(i)
				totalAngle = peakOrigin[2] + sampAngle
				xTotal = peakOrigin[0] - self.probe.segLength*cos(totalAngle)
				zTotal = peakOrigin[1] - self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				peakOrigin = pnt
				#peakCurves[-1].append(copy(peakOrigin))
				
				" if current joint is the maximum of a current peak, reset the origin "
				peakIndex = -1
				for index in range(len(self.localFit.jointClasses)):
					if self.localFit.jointClasses[index].count(i) > 0:
						peakIndex = index
						break
				
				if peakIndex >= 1 and max(self.localFit.jointClasses[peakIndex]) == i:
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
									
					
		errors = []
		for i in range(0,self.probe.numSegs-2):
			errors.append(fabs(self.probe.getServo(i)-self.probe.getServoCmd(i)))

		xP = []
		yP = []
		for p in actSegPoints:
			xP.append(p[0])
			yP.append(p[1])		
		pylab.plot(xP,yP,color='r')

		xP = []
		yP = []
		for p in self.cmdSegPoints:
			xP.append(p[0])
			yP.append(p[1])		
		pylab.scatter(xP,yP,color='g')
		
		for curve in peakCurves:
			xP = []
			yP = []
			for p in curve:
				xP.append(p[0])
				yP.append(p[1])
				
				
			pylab.scatter(xP,yP, linewidth=1, color='b')
			pylab.plot(xP,yP,color='b')
		
		for i in range(1, len(actSegPoints)-2):
			xP = [actSegPoints[i][0]]
			yP = [actSegPoints[i][1]]
			val = errors[i-1]
			if val > 1.0:
				val = 1.0
				
			pylab.scatter(xP,yP,color=str(1.0-val))

		pylab.xlim(-0.5,4)
		pylab.ylim(-2,2)
		pylab.savefig("fitPlot%04u.png" % self.plotCount)
		self.plotCount += 1
						

	def step(self):
		
		self.count += 1
		isDone = False

		" mask is used to determine the segments that need to be fitted to the curve "
		#self.mask = mask = [1.0 for i in range(39)]
		
		if self.transDone:
		#if self.count % 1000 == 0:
		
			self.transDone = False
			
			" compute the maximum error of the joints in the current peak "
			#peakErrorJoints = self.localFit.getPeakErrorJoints(self.currPeak)
			peakJoints = self.localFit.getPeakJoints(self.currPeak) + self.localFit.getPeakJoints(self.currPeak+1)
			peakJoints.sort()
			
			#print "peakJoints:", peakJoints
			
			errors = []
			for j in peakJoints:
				errors.append(self.probe.getServo(j)-self.probe.getServoCmd(j))
			
			#print "errors:", errors
			
			maxError = 0.0	
			for err in errors:
				if fabs(err) > maxError:
					maxError = fabs(err)
			
			" check for terminating criteria "
			currAmp = self.curve.getPeakAmp(self.currPeak)
			
			#currErr = self.localFit.getPeakError(self.currPeak)
			
			#print self.currPeak, "peakJoints =", peakJoints, "error =", maxError, currAmp, self.minAmp, self.maxAmp
			
			
			if len(peakJoints) == 0 and self.currPeak != 0:
				isDone = True
				print "terminating caseA"
				
			elif self.direction and 38 in peakJoints and currAmp > 0.5:
				isDone = True
				print "terminating caseB"
				
			elif not self.direction and 0 in peakJoints and currAmp > 0.5:
				isDone = True
				print "terminating caseC"
	
			if not isDone:
				
				" draw the state of the curve fitting "
				#self.drawFit()
				self.computeCmdSegPoints()

				currWidth = self.curve.getPeakWidth(self.currPeak)
				currAmp = self.curve.getPeakAmp(self.currPeak)
				nextVal = currAmp
				nextWidth = currWidth
	
					
				if not self.isJerking:
					" perform amplitude operation here "
					
					" error threshold that determines we have made an anchor contact "
					if maxError > 0.3 and self.maxAmp != 0.0:
						
						" if the difference between amplitudes is < 0.01 "
						if fabs(self.maxAmp-self.minAmp) < 0.01:
							
							" A: Section of snake successfully anchored "
							
							" set our new peak to solid so we don't recompute its values "
							if self.currPeak != 0:
								self.localFit.setSolid(self.currPeak)
								self.localFit.setSolid(self.currPeak+1)
								self.currPeak += 2
								self.minAmp = 0.0
								self.maxAmp = 0.0
								
								nextVal = 0.0
								self.ampInc = 0.04

								self.maxWidth = self.nomWidth
								self.minWidth = self.nomWidth

								nextWidth = self.nomWidth
								self.widthInc = 0.04
								self.isWidthMode = False
								self.widthLeft = False		
			
							else:
								
								self.isJerking = True
								

						else:
						
							" B: Section of snake experiencing error, now narrow down a tight fit"
						
							currAmp = self.curve.getPeakAmp(self.currPeak)
							
							if currAmp == self.minAmp:
								
								" C:  Overshoot.  Minimum is too high." 
								" Reduce amplitude minimum, set current to maximum"
								
								self.minAmp -= self.ampInc
								self.maxAmp = currAmp
								nextVal = self.minAmp
								
								#if currAmp >= 0.36:
								#	self.maxWidth = currWidth
								#	currWidth = self.minWidth
							
							else:
								
								"D: Bring amplitude down to minimum, and lower the step size"
								
								" bring it down to the minimum again "
								nextVal = self.minAmp
	
								#if currAmp >= 0.36:
								#	currWidth = self.minWidth
								
								" cut down the step size "
								self.ampInc /= 2.0
						
	
					else:
						
						"E:  Not experience any error, increase the amplitude by step size "

						" error not maxed out, so set min to the current amplitude "
						self.minAmp = currAmp
	
						nextVal = self.minAmp + self.ampInc
						
						" maximum is the value we just set "
						self.maxAmp = nextVal
			
						self.widthLeft = False		
						
				if self.isJerking:
					
					peakJoints = self.localFit.getPeakJoints(self.currPeak+1)
					if len(peakJoints) > 0:
						
						if self.direction:
							maxJoint = max(peakJoints)
							self.jerkJoint = maxJoint + 1
						else:
							minJoint = min(peakJoints)
							self.jerkJoint = minJoint - 1
						
						if self.jerkAngle == 0 and self.prevJerkAngle == 0:
							self.nomJerk = self.probe.getServo(self.jerkJoint) * 180.0 / pi
						
							print "nominal = " , self.nomJerk, self.probe.getServo(self.jerkJoint)
							self.prevJerkAngle = self.jerkAngle
							self.jerkAngle = 60	
						
						elif self.jerkAngle == 60:
							
							#print "setting jerk joint to ", self.jerkAngle + self.nomJerk
							print "jerk angle error = " , self.probe.getServo(j)-self.probe.getServoCmd(j)
							self.jerkErrors.append(self.probe.getServo(j)-self.probe.getServoCmd(j))
							self.prevJerkAngle = self.jerkAngle
							self.jerkAngle = -60
						
						elif self.jerkAngle == -60:
							print "jerk angle error = " , self.probe.getServo(j)-self.probe.getServoCmd(j)
							self.jerkErrors.append(self.probe.getServo(j)-self.probe.getServoCmd(j))

							" error = (-0.09, 0.005) is good "
							" error = (-2.73, -1.99) is bad "
							self.prevJerkAngle = self.jerkAngle
							self.jerkAngle = 0
						
						else:
							self.prevJerkAngle = self.jerkAngle
							self.jerkingDone = True

						#print "setting jerk joint to ", self.jerkAngle + self.nomJerk
						#print "jerk angle error = " , self.probe.getServo(j)-self.probe.getServoCmd(j)
			
				if self.isJerking and self.jerkingDone:
					
					" if the anchor is not secure, lets try it again "
					if abs(self.jerkErrors[0]) < 0.01 and abs(self.jerkErrors[1]) < 0.01:
						
						" reset the amplitude of peaks 0 and 1, nextVal = 0.0 "
						
						" increment the head length of adaptive cosine curve "
						head_len = self.curve.getHeadLength()
						self.curve.setHeadLength(head_len + 0.5)
						
						#tail_len = self.curve.getTailLength()
						#self.curve.setTailLength(tail_len - 0.5)
						
					else:
							
						self.localFit.setSolid(0)
						self.localFit.setSolid(1)
						self.currPeak += 2

						" reset the head length of adaptive cosine curve to 0 "
						#self.curve.setHeadLength(0.0)
						
					self.minAmp = 0.0
					self.maxAmp = 0.0

					nextVal = 0.0
					self.ampInc = 0.04

					self.maxWidth = self.nomWidth
					self.minWidth = self.nomWidth

					nextWidth = self.nomWidth
					self.widthInc = 0.04
					self.isWidthMode = False
					self.widthLeft = False		
					
					self.isJerking = False
					self.jerkingDone = False
					
					self.jerkErrors = []


				" reset all torques of the joints to maximum "
				for i in range(self.probe.numSegs-1):
					self.probe.setJointTorque(i, self.probe.maxTorque)

				" increase the amplitude depending on if this is the first peak or not "
				if self.currPeak == 0:					
					self.curve.setPeakAmp(0,nextVal)
					self.curve.setPeakAmp(1,nextVal)

				else:
					" stretch tail to remove closed-chain interference "
					self.curve.setTailLength(5.0)
					self.curve.setPeakAmp(self.currPeak,nextVal)
					self.curve.setPeakAmp(self.currPeak+1,nextVal)

					" weaken the feeder joints "
					" 30*2.5 is maximum "
					#peakJoints = self.localFit.getPeakJoints(self.currPeak)
					peakJoints = self.localFit.getPeakJoints(self.currPeak+1)
					if len(peakJoints) > 0:
						
						" weaken the joints so we don't move the snake's anchor's while stretching "
						if self.direction:
							maxJoint = max(peakJoints)
							for i in range(maxJoint+1, self.probe.numSegs-2):
								if i <= self.probe.numSegs-2:
									self.probe.setJointTorque(i, 3.0)
						else:
							minJoint = min(peakJoints)
							for i in range(0, minJoint):					
								if i <= self.probe.numSegs-2:
									self.probe.setJointTorque(i, 3.0)
									
				" execute the local curve fitting "
				self.localFit.step()
				self.holdT.reset(self.localFit.getJoints())
				
				if self.isJerking:
					print "setting jerk joint to ", self.nomJerk + self.jerkAngle
					self.holdT.positions[self.jerkJoint] = self.nomJerk + self.jerkAngle
					
				
				self.refDone = False
				self.transDone = self.holdT.step()
			else:
				# isDone = True
				
				" prevent the back anchor from slithering around or buckling "
				if self.direction:
					for i in range(self.probe.numSegs-8, self.probe.numSegs-2):				
						self.holdT.positions[i] = 180.0 / pi * self.probe.getServo(i)
				else:
					for i in range(0, 6):					
						self.holdT.positions[i] = 180.0 / pi * self.probe.getServo(i)

				self.holdT.step()
				
		else:
			
			" make no movements until all reference nodes are activated "
			if self.refDone:
				self.transDone = self.holdT.step()
	
		" collect joint settings for both behaviors "
		joints2 = self.localFit.getJoints()
		solidJoints = self.localFit.solidCompJoints
		
		joints3 = self.holdT.getJoints()


		#print joints3
		" Set the joints, with joints2 as priority over joints "
		if self.refDone:
			self.mergeJoints([joints3])
			
			" update mask for ContactReference "
			for i in range(self.probe.numSegs-1):
				if joints2[i] == None or self.localFit.isJointSolid(i):
					self.mask[i] = 1.0
				else:
					self.mask[i] = 0.0
		else:
			" transition the masks, and wait for them to be stable, then change the joints "

			allActive = True
			" update mask for ContactReference "
			for i in range(self.probe.numSegs-1):
				if joints2[i] == None or self.localFit.isJointSolid(i):
					self.mask[i] = 1.0

					" check if all reference points have been created "
					if not self.contacts.activeRef[i]:
						allActive = False
					
			if allActive:
				self.refDone = True
			
		
		
		if isDone:
			
			self.computeCenterPoints()
			self.curve.saveAmps()
			
			"reset everything "
			
			self.localFit = 0
			self.curve = 0
			self.computeCurve()
			
			self.mask = [0.0 for i in range(0,40)]
			self.count = 0
		
			self.currPeak = 0
		
		return isDone
	

	def reset(self):

		#self.computeCenterPoints()
		#self.curve.saveAmps()
		
		"reset everything "
		
		self.localFit = 0
		self.curve = 0
		self.computeCurve()
		
		self.mask = [0.0 for i in range(0,40)]
		self.count = 0
	
		self.currPeak = 0
	
	
	
	
