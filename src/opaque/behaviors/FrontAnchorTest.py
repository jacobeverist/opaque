import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from copy import *
from common import *
from Behavior import *
from FrontAnchorFit import FrontAnchorFit
from HoldTransition import HoldTransition
from AdaptiveAnchorCurve import AdaptiveAnchorCurve


class FrontAnchorTest(Behavior):

	def __init__(self, probe, contacts, mapGraph, direction = True):
		Behavior.__init__(self, probe)

		self.mapGraph = mapGraph
		self.localNode = self.mapGraph.getCurrentNode()
		self.contacts = contacts
		self.frontAnchorFit = 0
		self.concertinaFit = 0

		self.direction = direction
		self.isInit = False

		self.holdT = HoldTransition(probe)

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
			self.frontAnchorFit.setSide(True)
			
			self.concertinaFit.setSide(True)
			self.concertinaFit.setRootNode(38)
		else:
			self.frontAnchorFit.setSide(False)
			
			self.concertinaFit.setSide(False)
			self.concertinaFit.setRootNode(0)
		
	def computeCurve(self):
		
		self.frontCurve = AdaptiveAnchorCurve(4*pi)				

		if self.frontAnchorFit == 0:
			self.frontAnchorFit = FrontAnchorFit(self.probe, self.direction, 10)
			self.frontAnchorFit.setCurve(self.frontCurve)
		else:
			self.frontAnchorFit.setCurve(self.frontCurve)

		self.adaptiveCurve = AdaptiveCosine(4*pi)

		if self.concertinaFit == 0:
			self.concertinaFit = FastLocalCurveFit(self.probe, self.direction, 10)
			self.concertinaFit.setCurve(self.adaptiveCurve)
		else:
			self.concertinaFit.setCurve(self.adaptiveCurve)


		self.holdT.reset(self.frontAnchorFit.getJoints())
		self.refDone = False
			
	def getCenterPoints(self):
		
		#return copy(self.cmdSegPoints)
		return self.cPoints
			
	def computeCenterPoints(self):
		
		points = [[self.adaptiveCurve.ampSpacing/2.0 + i * self.adaptiveCurve.ampSpacing, -self.frontCurve.peakAmp[0]] for i in range(self.frontCurve.numPeaks)]

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
			
			#pose = self.contacts.getClosestPose(joint)
			pose = self.localNode.getJointPose(joint)
			point = [pose[0] + segLen*cos(pose[2]), pose[1] + segLen*sin(pose[2])]
			
			if point not in self.cPoints:
				self.cPoints.append(point)

	
	def computeCmdSegPoints(self):	
		
		if self.direction:
				
			" the configuration of the snake as we've commanded it "
			#self.cmdSegPoints = [[-self.probe.segLength,0.0,0.0]]
			self.cmdSegPoints = []
			cmdOrigin = [0.0,0.0,pi]
			self.cmdSegPoints.append(cmdOrigin)
			
			
			ind = range(0,10+1)
			ind.reverse()
			
			for i in ind:		
				sampAngle = self.probe.getServoCmd(i)
				totalAngle = cmdOrigin[2] + sampAngle
				xTotal = cmdOrigin[0] - self.probe.segLength*cos(totalAngle)
				zTotal = cmdOrigin[1] - self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				cmdOrigin = pnt
				self.cmdSegPoints.insert(0, cmdOrigin)

			cmdOrigin = [0.0,0.0,pi]
			for i in range(10,self.probe.numSegs-1):
				
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
		
		#peakJoints = self.frontAnchorFit.getPeakJoints(self.currPeak)

		pylab.clf()
		crvPoints = self.frontCurve.getPoints()		
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
			actSegPoints = []
			actOrigin = [0.0,0.0,pi]
			actSegPoints.append(actOrigin)

			ind = range(0,10+1)
			ind.reverse()
			
			for i in ind:
				sampAngle = self.probe.getServo(i)
				totalAngle = actOrigin[2] + sampAngle
				xTotal = actOrigin[0] - self.probe.segLength*cos(totalAngle)
				zTotal = actOrigin[1] - self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				actOrigin = pnt
				actSegPoints.insert(0,actOrigin)

			actOrigin = [0.0,0.0,pi]
			for i in range(11,self.probe.numSegs-1):
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

		"""
		for i in range(1, len(actSegPoints)-2):
			xP = [actSegPoints[i][0]]
			yP = [actSegPoints[i][1]]
			val = errors[i-1]
			if val > 1.0:
				val = 1.0
				
			pylab.scatter(xP,yP,color=str(1.0-val))
		"""

		pylab.xlim(-4,4)
		#pylab.xlim(-0.5,4)
		pylab.ylim(-3,3)
		pylab.savefig("fitPlot%04u.png" % self.plotCount)
		self.plotCount += 1
						

	def step(self):
		
		self.count += 1
		isDone = False

		" mask is used to determine the segments that need to be fitted to the curve "
		#self.mask = mask = [1.0 for i in range(39)]
		
		if self.transDone:
		
			self.transDone = False
			
			" compute the maximum error of the joints in the current peak "
			peakJoints = self.frontAnchorFit.getPeakJoints()
			peakJoints.sort()
			
			errors = []
			for j in peakJoints:
				errors.append(self.probe.getServo(j)-self.probe.getServoCmd(j))

			torques = []
			maxPeak = -1
			if len(peakJoints) > 0:
				maxPeak = max(peakJoints)

			for j in range(0,maxPeak+1):
				torques.append(self.probe.getJointTorque(j))
						
			maxError = 0.0	
			for err in errors:
				if fabs(err) > maxError:
					maxError = fabs(err)
			
			" check for terminating criteria "
			currAmp = self.frontCurve.getPeakAmp()
						
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
				self.drawFit()
				self.computeCmdSegPoints()

				currAmp = self.frontCurve.getPeakAmp()
				nextVal = currAmp	
				
				if not self.isJerking:
					" perform amplitude operation here "
					
					" error threshold that determines we have made an anchor contact "
					if maxError > 0.3 and self.maxAmp != 0.0:
						
						print errors
						print torques
						
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
					
				if self.isJerking:
					
					peakJoints = self.frontAnchorFit.getPeakJoints() 
					if len(peakJoints) > 0:
						
						if self.direction:
							maxJoint = max(peakJoints)
							self.jerkJoint = maxJoint + 1
							
							self.jerkJoints = []
							for k in range(0,4):
								jointNum = maxJoint + k + 1
								if jointNum > 38:
									break
								self.jerkJoints.append(jointNum)
								
						else:
							minJoint = min(peakJoints)
							self.jerkJoint = minJoint - 1
							self.jerkJoints = []
							for k in range(0,4):
								jointNum = minJoint - k - 1
								if jointNum < 0:
									break
								self.jerkJoints.append(jointNum)
													
						if self.jerkAngle == 0 and self.prevJerkAngle == 0:
							self.nomJerk = self.probe.getServo(self.jerkJoint) * 180.0 / pi
						
							for k in range(len(self.jerkJoints)):
								self.nomJerks[k] = self.probe.getServo(self.jerkJoints[k]) * 180.0 / pi
						
							print "nominal = " , self.nomJerk, self.probe.getServo(self.jerkJoint)
							self.prevJerkAngle = self.jerkAngle
							self.jerkAngle = 60	
						
							for k in range(len(self.jerkJoints)):
								self.prevJerkAngles[k] = self.jerkAngles[k]
								self.jerkAngles[k] = 60
								
							self.originPoses = []
							for k in peakJoints:
								self.originPoses.append(self.contacts.getAveragePose(k))
	
						elif self.jerkAngle == 60:
							
							#print "setting jerk joint to ", self.jerkAngle + self.nomJerk
							print "jerk angle error = " , self.probe.getServo(self.jerkJoint)-self.probe.getServoCmd(self.jerkJoint)

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
							for k in range(len(peakJoints)):
								tempPose = self.contacts.getAveragePose(peakJoints[k])
								originPose = self.originPoses[k]
								self.errorPoses1.append([tempPose[0]-originPose[0],tempPose[1]-originPose[1],tempPose[2]-originPose[2]])
						
							print self.errorPoses1
							
							self.originPoses = []
							for k in peakJoints:
								self.originPoses.append(self.contacts.getAveragePose(k))
							
						elif self.jerkAngle == -60:
							print "jerk angle error = " , self.probe.getServo(self.jerkJoint)-self.probe.getServoCmd(self.jerkJoint)
							self.jerkErrors.append(self.probe.getServo(self.jerkJoint)-self.probe.getServoCmd(self.jerkJoint))

							" error = (-0.09, 0.005) is good "
							" error = (-2.73, -1.99) is bad "
							self.prevJerkAngle = self.jerkAngle
							self.jerkAngle = 0
							
							for k in range(len(self.jerkJoints)):
								self.prevJerkAngles[k] = self.jerkAngles[k]
								self.jerkAngles[k] = 0

							self.errorPoses2 = []
							for k in range(len(peakJoints)):
								tempPose = self.contacts.getAveragePose(peakJoints[k])
								originPose = self.originPoses[k]
								self.errorPoses2.append([tempPose[0]-originPose[0],tempPose[1]-originPose[1],tempPose[2]-originPose[2]])
						
							print self.errorPoses2
						
						else:
							self.prevJerkAngle = self.jerkAngle
							self.jerkingDone = True

							for k in range(len(self.jerkJoints)):
								self.prevJerkAngles[k] = self.jerkAngles[k]

						#print "setting jerk joint to ", self.jerkAngle + self.nomJerk
						#print "jerk angle error = " , self.probe.getServo(j)-self.probe.getServoCmd(j)
			
				if self.isJerking and self.jerkingDone:
					
					" if the anchor is not secure, lets try it again "
					" 0.03707 "
					
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
					
					print "anchor errors =", err1, err2
					
					
					#if abs(self.jerkErrors[0]) < 0.01 and abs(self.jerkErrors[1]) < 0.01:
					#if abs(self.jerkErrors[0]) < 0.1 or abs(self.jerkErrors[1]) < 0.1:
					if err1 > 0.1 or err2 > 0.1:
					
						" reset the amplitude of peaks 0 and 1, nextVal = 0.0 "
						
						" increment the head length of adaptive cosine curve "
						head_len = self.frontCurve.getHeadLength()
						self.frontCurve.setHeadLength(head_len + 0.5)

						
						#tail_len = self.frontCurve.getTailLength()
						#self.frontCurve.setTailLength(tail_len - 0.5)
						
					else:
							
						self.currPeak += 2

						" reset the head length of adaptive cosine curve to 0 "
						#self.frontCurve.setHeadLength(0.0)
						
					self.minAmp = 0.0
					self.maxAmp = 0.0

					nextVal = 0.0
					self.ampInc = 0.04

	
					#nextVal = self.minAmp + self.ampInc
					" maximum is the value we just set "
					#self.maxAmp = nextVal
				
					self.isJerking = False
					self.jerkingDone = False
					
					self.jerkErrors = []


				" reset all torques of the joints to maximum "
				for i in range(self.probe.numSegs-1):
					self.probe.setJointTorque(i, self.probe.maxTorque)

				" increase the amplitude depending on if this is the first peak or not "
				if self.currPeak == 0:
					print "setting amp = " , nextVal		
					self.frontCurve.setPeakAmp(nextVal)

					" weaken head joints "
					" 30*2.5 is maximum "
					peakJoints = self.frontAnchorFit.getPeakJoints()

					print "peakJoints =", len(peakJoints), "amp =", self.frontCurve.getPeakAmp()
					if len(peakJoints) > 0 and self.frontCurve.getPeakAmp() > 0 and self.isJerking:
						
						print "direction =", self.direction
						if not self.direction:
							maxJoint = max(peakJoints)
							print "maxJoint = ", maxJoint
							for i in range(maxJoint+1, self.probe.numSegs-2):
								if i <= self.probe.numSegs-2:
									print "setting joint ", i, "to torque 3.0"
									self.probe.setJointTorque(i, 3.0)
						else:
							minJoint = min(peakJoints)
							print "minJoint = ", minJoint
							for i in range(0, minJoint):					
								if i <= self.probe.numSegs-2:
									print "setting joint ", i, "to torque 3.0"
									self.probe.setJointTorque(i, 3.0)

				else:
					" stretch tail to remove closed-chain interference "
					self.frontCurve.setTailLength(5.0)
					self.frontCurve.setPeakAmp(self.currPeak,nextVal)
					self.frontCurve.setPeakAmp(self.currPeak+1,nextVal)

					" weaken the feeder joints "
					" 30*2.5 is maximum "
					#peakJoints = self.frontAnchorFit.getPeakJoints(self.currPeak)
					peakJoints = self.frontAnchorFit.getPeakJoints(self.currPeak+1)
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
				self.frontAnchorFit.step()
				print "joints =", self.frontAnchorFit.getJoints()

				self.holdT.reset(self.frontAnchorFit.getJoints())
				
				if self.isJerking:
					print "setting jerk joint to ", self.nomJerk + self.jerkAngle
					#self.holdT.positions[self.jerkJoint] = self.nomJerk + self.jerkAngle
					for k in range(len(self.jerkJoints)):
						self.holdT.positions[self.jerkJoints[k]] = self.nomJerks[k] + self.jerkAngles[k]
					
				
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
		joints2 = self.frontAnchorFit.getJoints()

		joints3 = self.holdT.getJoints()


		#print joints2

		" Set the joints, with joints2 as priority over joints "
		if self.refDone:
			#print "joints =", joints3
			self.mergeJoints([joints3])
			
			" update mask for ContactReference "
			for i in range(self.probe.numSegs-1):
				if joints2[i] == None:
					self.mask[i] = 1.0
				else:
					self.mask[i] = 0.0
		else:
			" transition the masks, and wait for them to be stable, then change the joints "

			allActive = True
			" update mask for ContactReference "
			for i in range(self.probe.numSegs-1):
				if joints2[i] == None:
					self.mask[i] = 1.0

					" check if all reference points have been created "
					if not self.contacts.activeRef[i]:
						allActive = False
					
			if allActive:
				self.refDone = True
			
		
		
		if isDone:
			
			self.computeCenterPoints()
			
			"reset everything "
			
			self.frontAnchorFit = 0
			self.frontCurve = 0
			self.computeCurve()
			
			self.mask = [0.0 for i in range(0,40)]
			self.count = 0
		
			self.currPeak = 0
		
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
	
	
	
	
