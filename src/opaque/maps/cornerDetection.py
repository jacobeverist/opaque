import cv
from math import pi, sqrt, acos, cos, sin, fabs, floor
from numpy import arange, array
from copy import deepcopy, copy
from Pose import Pose

import Image

from scipy.cluster.vq import vq, kmeans, kmeans2, whiten
from scipy.cluster.hierarchy import fclusterdata
import pylab

saveCount = 0
histCount = 0

def Intersect(iSeg1, iSeg2):


	seg1 = [[float(iSeg1[0][0]), float(iSeg1[0][1])], [float(iSeg1[1][0]), float(iSeg1[1][1])]]
	seg2 = [[float(iSeg2[0][0]), float(iSeg2[0][1])], [float(iSeg2[1][0]), float(iSeg2[1][1])]]

	denom = (seg2[1][1] - seg2[0][1])*(seg1[1][0] - seg1[0][0]) - (seg2[1][0] - seg2[0][0])*(seg1[1][1] - seg1[0][1])
	nume_a = (seg2[1][0] - seg2[0][0])*(seg1[0][1] - seg2[0][1]) - (seg2[1][1] - seg2[0][1])*(seg1[0][0] - seg2[0][0])
	nume_b = (seg1[1][0] - seg1[0][0])*(seg1[0][1] - seg2[0][1]) - (seg1[1][1] - seg1[0][1])*(seg1[0][0] - seg2[0][0])
	intersection = [0.0,0.0]
	if denom == 0.0:
		if nume_a == 0.0 and nume_b == 0.0:
			return False, intersection
		return False, intersection

	ua = nume_a / denom
	ub = nume_b / denom

	if ua >= 0.0 and ua <= 1.0 and ub >= 0.0 and ub <= 1.0:
		# Get the intersection point.
		intersection[0] = seg1[0][0] + ua*(seg1[1][0] - seg1[0][0])
		intersection[1] = seg1[0][1] + ua*(seg1[1][1] - seg1[0][1])
		return True, (intersection[0],intersection[1])

	return False, (intersection[0],intersection[1])


def intersect2Lines(line1, line2):

	rho1 = line1[0]
	theta1 = line1[1]
	rho2 = line2[0]
	theta2 = line2[1]


	if theta1 == theta2:
		raise

	C1 = (rho2*sin(theta2) + ((rho2*cos(theta2) - rho1*cos(theta1))*cos(theta2))/sin(theta2) - rho1*sin(theta1)) / ( cos(theta1) - sin(theta1) * cos(theta2) / sin(theta2) )

	C2 = (rho2*cos(theta2) - rho1*cos(theta1) + C1*sin(theta1)) / sin(theta2)


	x = rho1*cos(theta1) - C1 * sin(theta1)
	y = rho1*sin(theta1) + C1 * cos(theta1)

	return (x, y)

def filterIntersections(lines, maxX, maxY):

	intersections = []

	BORDER = 0.2
	ANGLE_THRESH = pi / 8.

	limitX = BORDER*maxX
	limitY = BORDER*maxY

	for i in range(len(lines)):
		line1 = lines[i]
		for j in range(i+1, len(lines)):
			line2 = lines[j]

			theta1 = line1[1]
			theta2 = line2[1]

			try:
				pnt = intersect2Lines(lines[i], lines[j])

				if pnt[0] >= limitX and pnt[0] < maxX-limitX and pnt[1] >= limitY and pnt[1] < maxY-limitY:
					if theta1 > theta2:
						cornerAngle = theta1 - theta2
					else:
						cornerAngle = theta2 - theta1

					if cornerAngle > ANGLE_THRESH and cornerAngle < (pi - ANGLE_THRESH):
						intersections.append([pnt,line1,line2])
					#intersections.append([pnt,line1,line2])
			except:
				pass

	return intersections

def extractCornerCandidates(pilImg, estPose = []):

	global saveCount
	global histCount

	isRender = True
	



	" convert from PIL to OpenCV "
	#pilImg = Image.open('foo.png')       # PIL image
	img = cv.CreateImageHeader(pilImg.size, cv.IPL_DEPTH_8U, 1)
	cv.SetData(img, pilImg.tostring())

	#if isRender:
	#	cv.SaveImage("%04u_0.png" % saveCount, img)
	
	if isRender and len(estPose) > 0:
		color_dst = cv.CreateImage(cv.GetSize(img), 8, 3)
		cv.CvtColor(img, color_dst, cv.CV_GRAY2BGR)

		WLEN = 3.0
		WLEN2 = 5.0
		wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*cos(pi/3), -0.2 - WLEN*sin(pi/3)]]
		wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
		w1 = wall1[2]
		w2 = wall2[0]
		
		wall3 = [[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)], [0.4*cos(pi/6) - 4, 0.0], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
		lp = wall3[0]
		rp = wall3[2]
		
		wall6 = [lp, [lp[0] + WLEN*cos(pi/6), lp[1] + WLEN*sin(pi/6)]]
		wall6.append([wall6[1][0] + 0.4*cos(pi/3), wall6[1][1] - 0.4*sin(pi/3)])
		wall6.append([wall6[2][0] - WLEN*cos(pi/6), wall6[2][1] - WLEN*sin(pi/6)])
		wall6.append([wall6[3][0] + WLEN*cos(pi/3), wall6[3][1] - WLEN*sin(pi/3)])
		wall6.append([wall6[4][0] - 0.4*cos(pi/6), wall6[4][1] - 0.4*sin(pi/6)])
		wall6.append(w1)
		wall6.reverse()
		
		walls = [wall1, wall2, wall3, wall6]	
		for wall in walls:
			for i in range(len(wall)):
				p = copy(wall[i])
				p[0] += 6.0
				wall[i] = p		

		currProfile = Pose(estPose)
		
		divPix = 40.0
		numPixel = 401

		localWalls = []
		for wall in walls:
			newWall = []
			for i in range(len(wall)):
				p = copy(wall[i])
				localP = currProfile.convertGlobalToLocal(p)
				#indexX = localP[0]*divPix + numPixel/2 + 1
				#indexY = localP[1]*divPix + numPixel/2 + 1
				indexX = int(floor(localP[0]*divPix)) + numPixel/2 + 1
				indexY = int(floor(localP[1]*divPix)) + numPixel/2 + 1
				newWall.append([indexX, indexY])
			localWalls.append(newWall)
			
		for wall in localWalls:
			for i in range(len(wall)-1):
				indX = wall[i][0]
				indY = wall[i][1]
				indX_2 = wall[i+1][0]
				indY_2 = wall[i+1][1]
				cv.Line(color_dst, (indX, indY), (indX_2, indY_2), cv.CV_RGB(0, 0, 0), 1, 8)
		
		cv.SaveImage("%04u_0.png" % saveCount, color_dst)
	
	CORNER_THRESH = 400

	size = cv.GetSize(img)
	maxX = size[0]
	maxY = size[1]


	" 1st, run Sobel corner detector to detect corners "
	corners = cv.CreateImage(cv.GetSize(img), cv.IPL_DEPTH_32F, 1)
	cv.PreCornerDetect(img, corners, apertureSize=9)

	" 1st, run Sobel corner detector to detect corners "
	corners2 = cv.CreateImage(cv.GetSize(img), cv.IPL_DEPTH_32F, 1)
	cv.PreCornerDetect(img, corners2, apertureSize=13)


	" 2nd, run Canny edge detector to get the boundary points "
	edges = cv.CreateImage(cv.GetSize(img), cv.IPL_DEPTH_8U, 1)
	cv.Canny(img, edges, 50, 200, aperture_size=3)

	if isRender:
		cv.SaveImage("%04u_1.png" % saveCount, corners)
	#cv.SaveImage("%04u_2.png" % saveCount, edges)

	" 3rd, pick the maximum points out of the Sobel corner detection "
	count = 0
	maxVal = CORNER_THRESH
	maxPoints = []
	for m in range(maxX):
		for n in range(maxY):
			if corners[n,m] > maxVal:
				maxPoints.append((m,n))
				count += 1


	if isRender:	
		color_dst = cv.CreateImage(cv.GetSize(edges), 8, 3)
		cv.CvtColor(edges, color_dst, cv.CV_GRAY2BGR)
		for point in maxPoints:
			cv.Circle(color_dst, point, 2, cv.CV_RGB(0,0,255), thickness=1, lineType=8, shift=0) 
		cv.SaveImage("%04u_2.png" % saveCount, color_dst)
		
	#print count, "points over", maxVal

	#ROI_SIZE = 10
	ROI_SIZE = 15
	DIST_RES = 6
	ANG_RES = 0.1* pi /180.
	THRESH = 10 
	MAXLINEGAP = 3
	MINLINELENGTH = 3


	" print max points "
	filteredPoints = []

	cornerCandidates = []
	" 4th, for each maximum point, compute a corner centered on the point "
	for pnt in maxPoints:
		cv.SetImageROI(edges, (pnt[0]-(ROI_SIZE-1)/2, pnt[1]-(ROI_SIZE-1)/2, ROI_SIZE, ROI_SIZE))

		color_dst = cv.CreateImage(cv.GetSize(edges), 8, 3)
		cv.CvtColor(edges, color_dst, cv.CV_GRAY2BGR)

		" 5th, create a series of hypothesized lines going through the centroid "
		samples = []
		for theta in arange(-pi/2.-0.1, pi/2.+0.1, 0.05):
			rho = (ROI_SIZE-1.)/2. * cos(theta) + (ROI_SIZE-1.)/2. * sin(theta)
			samples.append([0, (rho, theta)])

		" 6th, for each point in the neighborhood, vote for each fitting line "
		size = cv.GetSize(edges)
		localMaxX = size[0]
		localMaxY = size[1]
		for m in range(localMaxX):
			for n in range(localMaxY):
				if edges[n,m] > 0:
					for k in range(len(samples)):
						rho1 = samples[k][1][0]
						theta1 = samples[k][1][1]
						rho2 = m * cos(theta1) + n * sin(theta1)
						dist = fabs(rho2-rho1)
						if dist <= 1.5:
							samples[k][0] += 1

		" 7th, select the top 10 most voted lines "
		results = deepcopy(samples)
		results.sort()
		results.reverse()


		#print "top edge votes for point:", pnt
		TOP_NUM = 10
		initVals = []
		for k in range(TOP_NUM):
			rho, theta = results[k][1]
			#print results[k][0]
			for m in range(results[k][0]):
				initVals.append(theta)
		#print


		" 8th, perform k-means clustering to select the two most likely edges that would form a corner "
		features  = array(initVals)
		features = features.reshape(len(initVals),1)
		#book1 = array([[pi/2.],[-pi/2.]])
		#book2 = array([[pi/4.],[-pi/4.]])
		#codebook1, distortion1 = kmeans2(features,book1, iter=100, minit='matrix')
		#codebook2, distortion2 = kmeans2(features,book2, iter=100, minit='matrix')
		resultVector = fclusterdata(features,t=2, criterion='maxclust')
		#theta1 = codebook2[0,0]
		#theta2 = codebook2[1,0]
		#rho1 = (ROI_SIZE-1.)/2. * cos(theta1) + (ROI_SIZE-1.)/2. * sin(theta1)
		#rho2 = (ROI_SIZE-1.)/2. * cos(theta2) + (ROI_SIZE-1.)/2. * sin(theta2)
		#selections = [(rho1,theta1), (rho2,theta2)]

		pylab.clf()
		blueVals = []
		redVals = []
		blueAvg = 0.0
		redAvg = 0.0
		for k in range(len(initVals)):
			val = features[k,0]
			clustID = resultVector[k]
			if clustID == 1:
				blueVals.append(val)
				blueAvg += val
			if clustID == 2:
				redVals.append(val)
				redAvg += val


		if len(redVals) > 0 and len(blueVals) > 0:

			blueAvg /= len(blueVals)
			redAvg /= len(redVals)

			theta1 = blueAvg
			theta2 = redAvg
			rho1 = (ROI_SIZE-1.)/2. * cos(theta1) + (ROI_SIZE-1.)/2. * sin(theta1)
			rho2 = (ROI_SIZE-1.)/2. * cos(theta2) + (ROI_SIZE-1.)/2. * sin(theta2)
			selections = [(rho1,theta1), (rho2,theta2)]

			#if len(redVals) > 0:
			#	pylab.hist(redVals, bins=arange(-pi/2.-0.1, pi/2.+0.1, 0.05),color='r')
			#if len(blueVals) > 0:
			#	pylab.hist(blueVals, bins=arange(-pi/2.-0.1, pi/2.+0.1, 0.05),color='b')
			#pylab.savefig("hist_%04u.png" % histCount)
			#histCount += 1


			" 9th, check if we've successfully clustered two results by checking for default answers "
			#isSuccess = True
			#if theta1 == pi/4. or theta1 == -pi/4.:
			#	if samples[-19][0] < 9:
			#		isSuccess = False
			#		maxDist = 0
			#		maxK = 0
			#		for k in range(10):
			#			if fabs(results[k][1][1] - theta2) > maxDist:
			#				maxDist = fabs(results[k][1][1] - theta2)
			#				maxK = k
			#		theta1 = results[maxK][1][1]
			#		rho1 = results[maxK][1][0]

			#if theta2 == -pi/4. or theta2 == pi/4.:
			#	if samples[18][0] < 9:
			#		isSuccess = False
			#		maxDist = 0
			#		maxK = 0
			#		for k in range(10):
			#			if fabs(results[k][1][1] - theta1) > maxDist:
			#				maxDist = fabs(results[k][1][1] - theta1)
			#				maxK = k
			#		theta2 = results[maxK][1][1]
			#		rho2 = results[maxK][1][0]

			isSuccess = True

			" 10th, check if corner is sufficiently acute to be used as a feature "
			isCorner = False
			if theta1 > theta2:
				cornerAngle = theta1 - theta2
			else:
				cornerAngle = theta2 - theta1

			ANGLE_THRESH = pi / 8.
			if cornerAngle > ANGLE_THRESH and cornerAngle < (pi - ANGLE_THRESH):
				isCorner = True

				#print "cornerAngle =", cornerAngle
			#isCorner = True

			if isSuccess and isCorner:


				#top = []
				#for k in range(10):
				#	rho, theta = results[k][1]
				#	top.append(theta)
				#	a = cos(theta)
				#	b = sin(theta)
				#	x0 = a * rho 
				#	y0 = b * rho
				#	pt1 = (cv.Round(x0 + 400*(-b)), cv.Round(y0 + 400*(a)))
				#	pt2 = (cv.Round(x0 - 400*(-b)), cv.Round(y0 - 400*(a)))
				#	seg = (pt1,pt2)
				#	cv.Line(color_dst, seg[0], seg[1], cv.CV_RGB(0, 0, 255), 1, 8)

				#for (rho,theta) in selections:
				#	a = cos(theta)
				#	b = sin(theta)
				#	x0 = a * rho 
				#	y0 = b * rho
				#	pt1 = (cv.Round(x0 + 400*(-b)), cv.Round(y0 + 400*(a)))
				#	pt2 = (cv.Round(x0 - 400*(-b)), cv.Round(y0 - 400*(a)))
				#	seg = (pt1,pt2)
				#	cv.Line(color_dst, seg[0], seg[1], cv.CV_RGB(255, 0, 0), 1, 8)
		
				#cv.SaveImage("%04u_1.png" % saveCount, edges)
				#cv.SaveImage("%04u_2.png" % saveCount, color_dst)
				#saveCount += 1

				" 11th, convert to full image coordinates " 
				rho1 = selections[0][0]
				theta1 = selections[0][1]
				rho2 = selections[1][0]
				theta2 = selections[1][1]

				globalRho1 = pnt[0] * cos(theta1) + pnt[1] * sin(theta1)
				globalRho2 = pnt[0] * cos(theta2) + pnt[1] * sin(theta2)

				point = pnt

				a = cos(theta1)
				b = sin(theta1)
				x0 = a * globalRho1 
				y0 = b * globalRho1
				pt1 = (cv.Round(x0 + 400*(-b)), cv.Round(y0 + 400*(a)))
				pt2 = (cv.Round(x0 - 400*(-b)), cv.Round(y0 - 400*(a)))
				vec1_1 = [pt1[0]-point[0], pt1[1]-point[1]]
				vec1_2 = [pt2[0]-point[0], pt2[1]-point[1]]
				mag1_1 = sqrt(vec1_1[0]**2 + vec1_1[1]**2)
				mag1_2 = sqrt(vec1_2[0]**2 + vec1_2[1]**2)
				vec1_1[0] /= mag1_1
				vec1_1[1] /= mag1_1
				vec1_2[0] /= mag1_2
				vec1_2[1] /= mag1_2

				a = cos(theta2)
				b = sin(theta2)
				x0 = a * globalRho2 
				y0 = b * globalRho2
				pt1 = (cv.Round(x0 + 400*(-b)), cv.Round(y0 + 400*(a)))
				pt2 = (cv.Round(x0 - 400*(-b)), cv.Round(y0 - 400*(a)))
				vec2_1 = [pt1[0]-point[0], pt1[1]-point[1]]
				vec2_2 = [pt2[0]-point[0], pt2[1]-point[1]]
				mag2_1 = sqrt(vec2_1[0]**2 + vec2_1[1]**2)
				mag2_2 = sqrt(vec2_2[0]**2 + vec2_2[1]**2)
				vec2_1[0] /= mag2_1
				vec2_1[1] /= mag2_1
				vec2_2[0] /= mag2_2
				vec2_2[1] /= mag2_2

				cornerAngle = fabs(theta1 - theta2)
				#print point[0], point[1], cornerAngle, vec1_1, vec1_2, vec2_1, vec2_2
				#print vec1_1, vec1_2, vec2_1, vec2_2

				" vec2_1, vec1_1 "
				dirA = [vec2_1[0]+vec1_1[0], vec2_1[1]+vec1_1[1]]
				magA = sqrt(dirA[0]**2 + dirA[1]**2)
				#print "dirA =", dirA
				#print "magA =", magA
				dirA[0] /= magA
				dirA[1] /= magA
				indX = cv.Floor(point[0] + 5*dirA[0])
				indY = cv.Floor(point[1] + 5*dirA[1])
				valA = img[indY,indX]
				" vec1_1, vec2_2 "
				dirB = [vec1_1[0]+vec2_2[0], vec1_1[1]+vec2_2[1]]
				magB = sqrt(dirB[0]**2 + dirB[1]**2)
				#print "dirB =", dirB
				#print "magB =", magB
				dirB[0] /= magB
				dirB[1] /= magB
				indX = cv.Floor(point[0] + 5*dirB[0])
				indY = cv.Floor(point[1] + 5*dirB[1])
				valB = img[indY,indX]
				" vec2_2, vec1_2 "
				dirC = [vec2_2[0]+vec1_2[0], vec2_2[1]+vec1_2[1]]
				magC = sqrt(dirC[0]**2 + dirC[1]**2)
				#print "dirC =", dirC
				#print "magC =", magC
				dirC[0] /= magC
				dirC[1] /= magC
				indX = cv.Floor(point[0] + 5*dirC[0])
				indY = cv.Floor(point[1] + 5*dirC[1])
				valC = img[indY,indX]
				" vec1_2, vec2_1 "
				dirD = [vec1_2[0]+vec2_1[0], vec1_2[1]+vec2_1[1]]
				magD = sqrt(dirD[0]**2 + dirD[1]**2)
				#print "dirD =", dirD
				#print "magD =", magD
				dirD[0] /= magD
				dirD[1] /= magD
				indX = cv.Floor(point[0] + 5*dirD[0])
				indY = cv.Floor(point[1] + 5*dirD[1])
				valD = img[indY,indX]
				
				statA = 0
				statB = 0
				statC = 0
				statD = 0
				if valA > 128:
					statA = 1
				if valB > 128:
					statB = 1
				if valC > 128:
					statC = 1
				if valD > 128:
					statD = 1


				#print valA, valB, valC, valD
				#print statA, statB, statC, statD

				if statA + statB + statC + statD == 3:

					cornerAngle = 0.0
					if statA == 0:
						" vec2_1, vec1_1 "
						inwardVec = dirC 
						#dotProduct = vec2_1[0]*vec1_1[0] + vec2_1[1]*vec1_1[1]
						dotProduct1 = inwardVec[0]*vec1_1[0] + inwardVec[1]*vec1_1[1]
						dotProduct2 = vec2_1[0]*inwardVec[0] + vec2_1[1]*inwardVec[1]
						#cornerAngle = acos(dotProduct)
						cornerAngle = acos(dotProduct1) + acos(dotProduct2)
					if statB == 0:
						" vec1_1, vec2_2 "
						inwardVec = dirD 
						#dotProduct = vec1_1[0]*vec2_2[0] + vec1_1[1]*vec2_2[1]
						dotProduct1 = inwardVec[0]*vec2_2[0] + inwardVec[1]*vec2_2[1]
						dotProduct2 = vec1_1[0]*inwardVec[0] + vec1_1[1]*inwardVec[1]
						#cornerAngle = acos(dotProduct)
						cornerAngle = acos(dotProduct1) + acos(dotProduct2)
					if statC == 0:
						" vec2_2, vec1_2 "
						inwardVec = dirA 
						#dotProduct = vec2_2[0]*vec1_2[0] + vec2_2[1]*vec1_2[1]
						dotProduct1 = inwardVec[0]*vec1_2[0] + inwardVec[1]*vec1_2[1]
						dotProduct2 = vec2_2[0]*inwardVec[0] + vec2_2[1]*inwardVec[1]
						#cornerAngle = acos(dotProduct)
						cornerAngle = acos(dotProduct1) + acos(dotProduct2)
					if statD == 0:
						" vec1_2, vec2_1 "
						inwardVec = dirB 
						#dotProduct = vec1_2[0]*vec2_1[0] + vec1_2[1]*vec2_1[1]
						dotProduct1 = inwardVec[0]*vec2_1[0] + inwardVec[1]*vec2_1[1]
						dotProduct2 = vec1_2[0]*inwardVec[0] + vec1_2[1]*inwardVec[1]
						#cornerAngle = acos(dotProduct)
						cornerAngle = acos(dotProduct1) + acos(dotProduct2)

					#print "cornerAngle =", cornerAngle
					cornerCandidates.append(((globalRho1,theta1),(globalRho2,theta2), pnt, inwardVec, cornerAngle))

	points = []
	for ((rho1,theta1),(rho2,theta2), point, inwardVec, cornerAngle) in cornerCandidates:
		points.append(point)

	finalCandidates = []
	features = array(points)
	if len(points) > 1:
		features = features.reshape(len(points),2)
		#print features
		resultVector = fclusterdata(features,t=4, criterion='maxclust')
		#print resultVector

		numCandidates = resultVector.max()
		avgs = [[[0.0,0.0],0.0,[0.0, 0.0],0.0,0.0] for i in range(numCandidates)]
		nums = [0 for i in range(numCandidates)]

		for k in range(len(points)):
			valX = features[k,0]
			valY = features[k,1]
			cornerVal = cornerCandidates[k][4]
			inwardVec = cornerCandidates[k][3]

			theta1 = cornerCandidates[k][0][1]
			theta2 = cornerCandidates[k][1][1]

			clustID = resultVector[k]
			" average point "
			avgs[clustID-1][0][0] += valX
			avgs[clustID-1][0][1] += valY

			" average corner angle "
			avgs[clustID-1][1] += cornerVal

			" average inward vector "
			avgs[clustID-1][2][0] += inwardVec[0]
			avgs[clustID-1][2][1] += inwardVec[1]

			" average the thetas "
			avgs[clustID-1][3] += theta1
			avgs[clustID-1][4] += theta2
			
			nums[clustID-1] += 1

		for k in range(numCandidates):
			avgs[k][0][0] /= nums[k]
			avgs[k][0][1] /= nums[k]
			avgs[k][1] /= nums[k]

			vec = avgs[k][2]
			mag = sqrt(vec[0]**2 + vec[1]**2)
			vec[0] /= mag
			vec[1] /= mag
			avgs[k][2] = vec

			avgs[k][3] /= nums[k]
			avgs[k][4] /= nums[k]


			finalCandidates.append((avgs[k][0], avgs[k][1], avgs[k][2], avgs[k][3], avgs[k][4]))
	else:

		for ((rho1,theta1),(rho2,theta2), point, inwardVec, cornerAngle) in cornerCandidates:
			finalCandidates.append((point, cornerAngle, inwardVec, theta1, theta2))

	cv.ResetImageROI(edges)


	ROI_SIZE_2 = 21
	finalFinalCandidates = []

	for (point, cornerAngle, inwardVec, theta1, theta2) in finalCandidates:
		pnt = point
		#cv.SetImageROI(edges, (pnt[0]-(ROI_SIZE_2-1)/2, pnt[1]-(ROI_SIZE_2-1)/2, ROI_SIZE_2, ROI_SIZE_2))

		" 5th, create a series of hypothesized lines going through the centroid "
		rho1 = (ROI_SIZE_2-1.)/2. * cos(theta1) + (ROI_SIZE_2-1.)/2. * sin(theta1)
		rho2 = (ROI_SIZE_2-1.)/2. * cos(theta2) + (ROI_SIZE_2-1.)/2. * sin(theta2)

		halfAngle = cornerAngle / 2.0
		rightVec = [inwardVec[0]*cos(halfAngle) - inwardVec[1]*sin(halfAngle), inwardVec[0]*sin(halfAngle) + inwardVec[1]*cos(halfAngle)]
		leftVec = [inwardVec[0]*cos(halfAngle) + inwardVec[1]*sin(halfAngle), inwardVec[0]*sin(halfAngle) - inwardVec[1]*cos(halfAngle)]

		distCount1 = 0
		distCount2 = 0
		
		WID = 2
		
		for i in range(0,21):
			ptR = [pnt[0] + 0.5 + rightVec[0]*i, pnt[1] + 0.5 + rightVec[1]*i]
			
			
			ptR[0] = int(floor(ptR[0]))
			ptR[1] = int(floor(ptR[1]))
			hasEdge = False
			for m in range(ptR[0]-WID,ptR[0]+WID+1):
				for n in range(ptR[1]-WID,ptR[1]+WID+1):
					#print m, n, edges[m,n]
					if edges[n,m] > 0:
						hasEdge = True
			if hasEdge:
				distCount1 += 1
			else:
				break
			
		for i in range(0,21):
			ptL = [pnt[0] + 0.5 + leftVec[0]*i, pnt[1] + 0.5 + leftVec[1]*i]

			ptL[0] = int(floor(ptL[0]))
			ptL[1] = int(floor(ptL[1]))
			hasEdge = False
			for m in range(ptL[0]-WID,ptL[0]+WID+1):
				for n in range(ptL[1]-WID,ptL[1]+WID+1):
					if edges[n,m] > 0:
						hasEdge = True
			if hasEdge:
				distCount2 += 1
			else:
				break
		
		"""
		if ptR[0] > pnt[0]:
			rangeX = range(pnt[0],ptR[0]+1)
		else:
			rangeX = range(ptR[0],pnt[0]+1)

		if ptR[1] > pnt[1]:
			rangeY = range(pnt[1],ptR[1]+1)
		else:
			rangeY = range(ptR[1],pnt[1]+1)
			
		for m in rangeX:
			for n in rangeY:
				if edges[n,m] > 0:
					rhoN = m * cos(theta1) + n * sin(theta1)
					dist = fabs(rhoN-rho1)
					if dist <= 1.0:
						distCount1 += 1

		" 6th, for each point in the neighborhood, vote for each fitting line "
		distCount1 = 0
		distCount2 = 0
		size = cv.GetSize(edges)
		localMaxX = size[0]
		localMaxY = size[1]
		for m in range(localMaxX):
			for n in range(localMaxY):
				if edges[n,m] > 0:
					rhoN = m * cos(theta1) + n * sin(theta1)
					dist = fabs(rhoN-rho1)
					if dist <= 1.0:
						distCount1 += 1

					rhoN = m * cos(theta2) + n * sin(theta2)
					dist = fabs(rhoN-rho2)
					if dist <= 1.0:
						distCount2 += 1
		"""

		
		cornerVal = corners2[pnt[1],pnt[0]]

		if cornerVal >= 1000000:
			finalFinalCandidates.append((point, cornerAngle, inwardVec, theta1, theta2))

		print "corner candidate:", cornerVal, distCount1, distCount2, pnt, cornerAngle, inwardVec
		print

	cv.ResetImageROI(edges)

	if isRender:
		color_dst = cv.CreateImage(cv.GetSize(edges), 8, 3)
		cv.CvtColor(edges, color_dst, cv.CV_GRAY2BGR)
		angles = []
		for (point, cornerAngle, inwardVec, theta1, theta2) in finalFinalCandidates:
	
			angles.append(cornerAngle)
	
			indX = cv.Floor(point[0])
			indY = cv.Floor(point[1])
	
			pt2 = [point[0] + inwardVec[0]*8, point[1] + inwardVec[1]*8]
			indX_2 = cv.Floor(pt2[0])
			indY_2 = cv.Floor(pt2[1])
			cv.Line(color_dst, (indX, indY), (indX_2, indY_2), cv.CV_RGB(255, 0, 0), 2, 8)
	
	
		numCorners = len(finalFinalCandidates)
		font = cv.InitFont(cv.CV_FONT_HERSHEY_COMPLEX, 0.5, 0.5, 0.0, 1, cv.CV_AA)
		cv.PutText(color_dst, "%d candidate corners" % numCorners, (10, 50), font, cv.RGB(255,255,255)) 
		cv.SaveImage("%04u_3.png" % saveCount, color_dst)
		saveCount += 1

	" compute the direction of the corner.  vector goes inward towards open space "

	#print len(finalCandidates), "corner candidates"
	return finalFinalCandidates

def findCornerCandidates(lines, edgeImage, cornersImage):

	size = cv.GetSize(edgeImage)
	maxX = size[0]
	maxY = size[1]

	intersections = filterIntersections(lines, maxX, maxY)


	count = 0
	maxVal = 400
	for i in range(maxX):
		for j in range(maxY):
			if cornersImage[j,i] > maxVal:
				count += 1
	#print count, "points over", maxVal


	" remove non-corner candidates "
	candidates = []
	for inter in intersections:
		point = inter[0]
		indX = cv.Floor(point[0])
		indY = cv.Floor(point[1])

		if cornersImage[indY, indX] > 0.0:
			candidates.append(inter)

		

	#SAMPLE_COUNT = 20
	SAMPLE_COUNT = 7
	APERT = 1
	HIGH_THRESH = 6
	LOW_THRESH = 1

	results = []

	for inter in candidates:

		left1Count = 0
		right1Count = 0
		left2Count = 0
		right2Count = 0

		point = inter[0]
		line1 = inter[1]
		line2 = inter[2]

		rho1 = line1[0]
		theta1 = line1[1]

		rho2 = line2[0]
		theta2 = line2[1]

		#print "point:", point
		#point[0] = rho1*cos(theta1) - C1 * sin(theta1)
		#point[1] = rho1*sin(theta1) + C1 * cos(theta1)

		C1 = 0
		try:
			C1 = (point[1] - rho1*sin(theta1)) / cos(theta1)
		except:
			C1 =  -(point[0] - rho1*cos(theta1)) / sin(theta1)

		C2 = 0
		try:
			C2 = (point[1] - rho2*sin(theta2)) / cos(theta2)
		except:
			C2 = -(point[0] - rho2*cos(theta2)) / sin(theta2)

		#print "rho1, theta1, rho2, theta2 =", rho1, theta1, rho2, theta2
		#print "C1,C2 =", C1, C2

		" C1 is increasing "
		for i in range(SAMPLE_COUNT):
			sample = [0., 0.]
			sample[0] = rho1*cos(theta1) - (C1+i) * sin(theta1)
			sample[1] = rho1*sin(theta1) + (C1+i) * cos(theta1)

			indX = cv.Floor(sample[0])
			indY = cv.Floor(sample[1])

			isEdge = False
			for m in range(-APERT,APERT+1):
				for n in range(-APERT,APERT+1):
					if edgeImage[(indY+n,indX+m)] > 0:
						isEdge = True
			if isEdge:
				left1Count += 1

		" C1 is decreasing "
		for i in range(SAMPLE_COUNT):
			sample = [0., 0.]
			sample[0] = rho1*cos(theta1) - (C1-i) * sin(theta1)
			sample[1] = rho1*sin(theta1) + (C1-i) * cos(theta1)
			
			indX = cv.Floor(sample[0])
			indY = cv.Floor(sample[1])

			isEdge = False
			for m in range(-APERT,APERT+1):
				for n in range(-APERT,APERT+1):
					if edgeImage[(indY+n,indX+m)] > 0:
						isEdge = True
			if isEdge:
				right1Count += 1

		" C2 is increasing "
		for i in range(SAMPLE_COUNT):
			sample = [0., 0.]
			sample[0] = rho2*cos(theta2) - (C2+i) * sin(theta2)
			sample[1] = rho2*sin(theta2) + (C2+i) * cos(theta2)

			indX = cv.Floor(sample[0])
			indY = cv.Floor(sample[1])

			isEdge = False
			for m in range(-APERT,APERT+1):
				for n in range(-APERT,APERT+1):
					if edgeImage[(indY+n,indX+m)] > 0:
						isEdge = True

			if isEdge:
				left2Count += 1

		" C2 is decreasing "
		for i in range(SAMPLE_COUNT):
			sample = [0., 0.]
			sample[0] = rho2*cos(theta2) - (C2-i) * sin(theta2)
			sample[1] = rho2*sin(theta2) + (C2-i) * cos(theta2)

			indX = cv.Floor(sample[0])
			indY = cv.Floor(sample[1])

			for m in range(-APERT,APERT+1):
				for n in range(-APERT,APERT+1):
					if edgeImage[(indY+n,indX+m)] > 0:
						isEdge = True

			if isEdge:
				right2Count += 1
				
		#results.append(point)

		if (left1Count >= HIGH_THRESH and right1Count <= LOW_THRESH) or (left1Count <= LOW_THRESH and right1Count >= HIGH_THRESH):
			if (left2Count >= HIGH_THRESH and right2Count <= LOW_THRESH) or (left2Count <= LOW_THRESH and right2Count >= HIGH_THRESH):
				#print left1Count, right1Count, left2Count, right2Count
				results.append([point, line1, line2])

		#if (left1Count > 0 or right2Count > 0) and (left2Count > 0 or right2Count > 0):

	return results

if __name__ == '__main__':
	
	count = 0
	total = 120
	
	dirName = "../testData/sensorLocalize2"
	fileName = dirName + "/localSweepMap%03u_0000.png"			
	
	#positiveExamples = [19]
	positiveExamples = [2,4,17,19,27,35,37,39,44,45,51,53,56,59]
	#positiveExamples = [37,45,51,53]
	
	#for i in positiveExamples:
	for i in range(total):
		fileName % i
		img0 = cv.LoadImage(fileName % i, cv.CV_LOAD_IMAGE_UNCHANGED)
	
		extractCornerCandidates(img0)
	
	
		"""
	
		size = cv.GetSize(img0)
		size = (6*size[0], size[1])
		#eig_image = cv.CreateImage(size, cv.IPL_DEPTH_32F, 1)
		harris_dst = cv.CreateImage(cv.GetSize(img0), cv.IPL_DEPTH_32F, 1)
		corners = cv.CreateImage(cv.GetSize(img0), cv.IPL_DEPTH_32F, 1)
		edges = cv.CreateImage(cv.GetSize(img0), cv.IPL_DEPTH_8U, 1)
		#temp_image = cv.CreateImage(cv.GetSize(img0), cv.IPL_DEPTH_8U, 1)
	
		#cv.CornerEigenValsAndVecs(img0, eig_image, 30, aperture_size=9)
		#cv.CornerHarris(img0, harris_dst, 10, aperture_size=9, k=0.04)
		cv.CornerHarris(img0, harris_dst, 30, aperture_size=9, k=0.1)
		cv.Canny(img0, edges, 50, 200, aperture_size=3)
	
		edgeCopy = cv.CloneImage(edges)
	
		#lines = cv.HoughLines2(edges, cv.CreateMemStorage(), cv.CV_HOUGH_PROBABILISTIC, 5, pi/180., 50, 50, 10)
		#lines = cv.HoughLines2(edges, cv.CreateMemStorage(), cv.CV_HOUGH_PROBABILISTIC, 1, pi / 180, 10, 10, 1)
	
		DIST_RES = 1
		ANG_RES = 1* pi /180.
		MAXLINEGAP = 4
		#MAXLINEGAP = 10
		#MINLINELENGTH = 10
		MINLINELENGTH = 5
		#THRESH = 10
		THRESH = 10
		USE_STANDARD = True
		
		color_dst = cv.CreateImage(cv.GetSize(img0), 8, 3)
		#cv.SaveImage("%04u_4.png" % i, edges)
		cv.CvtColor(edges, color_dst, cv.CV_GRAY2BGR)
		newLines = []
		polarLines = []
		filteredLines = []
	
		if USE_STANDARD:
			storage = cv.CreateMemStorage()
			polarLines = cv.HoughLines2(edges, storage, cv.CV_HOUGH_STANDARD, DIST_RES, ANG_RES, THRESH, 0, 0)
	
	
			#print len(polarLines), "lines for image", i
			for (rho, theta) in polarLines[:100]:
				#print "theta:", theta
				a = cos(theta)
				b = sin(theta)
				x0 = a * rho 
				y0 = b * rho
				pt1 = (cv.Round(x0 + 400*(-b)), cv.Round(y0 + 400*(a)))
				pt2 = (cv.Round(x0 - 400*(-b)), cv.Round(y0 - 400*(a)))
				newLines.append((pt1,pt2))
				#cv.Line(color_dst, pt1, pt2, cv.RGB(255, 0, 0), 1, 8)
		else:
			newLines = cv.HoughLines2(edges, cv.CreateMemStorage(), cv.CV_HOUGH_PROBABILISTIC, DIST_RES, ANG_RES, THRESH, MINLINELENGTH, MAXLINEGAP)
			#cv.Line(color_dst, line[0], line[1], cv.CV_RGB(255, 0, 0), 2, 8)
	
		points = []
		size = cv.GetSize(edgeCopy)
		maxX = size[0]
		maxY = size[1]
	
		#print "size:", maxX, maxY
		for j in range(len(newLines)):
			line1 = newLines[j]
			for k in range(j+1, len(newLines)):
				line2 = newLines[k]
				p1_1 = line1[0]
				p1_2 = line1[1]
				p2_1 = line2[0]
				p2_2 = line2[1]
				result, pnt = Intersect(line1, line2)
	
	
				if result and pnt[0] >= 0 and pnt[0] < maxX and pnt[1] >= 0 and pnt[1] < maxY:
	
	
					try:
						indX = cv.Floor(pnt[0])
						indY = cv.Floor(pnt[1])
	
						isEdge = False
						#for m in range(-1,2):
						#	for n in range(-1,2):
						#		if edgeCopy[(indY+n,indX+m)] > 0:
						#			isEdge = True
						#			#print "edge at:", indX+m, indY+n, edgeCopy[indY+n,indX+m]
						if edgeCopy[(indY,indX)] > 0:
							isEdge = True
	
						if isEdge:
							vec1 = [p1_1[0]-p1_2[0], p1_1[1]-p1_2[1]]
							vec2 = [p2_1[0]-p2_2[0], p2_1[1]-p2_2[1]]
							dotProduct = vec1[0]*vec2[0] + vec1[1]*vec2[1]
							mag1 = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1])
							mag2 = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1])
							try:
								cornerAngle = acos(dotProduct/(mag1*mag2))
								#if cornerAngle > 0.25 and cornerAngle < 2.5:
								if cornerAngle > pi/4.0 and cornerAngle < 3.*pi/4.:
									points.append((indX,indY))
									filteredLines.append(line1)
									filteredLines.append(line2)
									#print line1, line2, points[-1], edgeCopy[(indY,indX)]
									#cv.Circle(color_dst, (indX,indY), 5, cv.CV_RGB(0,0,255), thickness=1, lineType=8, shift=0) 
							except:
								pass
	
	
						#if edges[cv.Floor(pnt[0]), cv.Floor(pnt[1])] > 0:
						#	print edges[cv.Floor(pnt[0]), cv.Floor(pnt[1])] 
					except:
						#print "FAILED:", pnt
						pass
	
	
	
				#if result and pnt[0] < size[0] and pnt[0] >= 0 and pnt[1] < size[1] and pnt[1] >= 0:
				#	print pnt
				#	points.append(pnt)
	
		"""
	
		"""
		EDGE_DIST = 10.0
		numLines = len(lines)
		toDraw = []
		points = []
		for j in range(len(lines)):
			line1 = lines[j]
			for k in range(j+1, len(lines)):
				line2 = lines[k]
				p1_1 = line1[0]
				p1_2 = line1[1]
				p2_1 = line2[0]
				p2_2 = line2[1]
	
				dist1 = sqrt((p1_1[0]-p2_1[0])**2 + (p1_1[1]-p2_1[1])**2)
				dist2 = sqrt((p1_2[0]-p2_1[0])**2 + (p1_2[1]-p2_1[1])**2)
				dist3 = sqrt((p1_2[0]-p2_2[0])**2 + (p1_2[1]-p2_2[1])**2)
	
	
				if dist1 <= EDGE_DIST or dist2 <= EDGE_DIST or dist3 <= EDGE_DIST:
					if dist1 <= dist2 and dist1 <= dist3:
						vec1 = [p1_2[0]-p1_1[0], p1_2[1]-p1_1[1]]
						vec2 = [p2_2[0]-p2_1[0], p2_2[1]-p2_1[1]]
						pnt = p1_1
	
					elif dist2 <= dist1 and dist2 <= dist3:
						vec1 = [p1_1[0]-p1_2[0], p1_1[1]-p1_2[1]]
						vec2 = [p2_2[0]-p2_1[0], p2_2[1]-p2_1[1]]
						pnt = p1_2
	
					elif dist3 <= dist1 and dist3 <= dist2:
						vec1 = [p1_1[0]-p1_2[0], p1_1[1]-p1_2[1]]
						vec2 = [p2_1[0]-p2_2[0], p2_1[1]-p2_2[1]]
						pnt = p1_2
	
					dotProduct = vec1[0]*vec2[0] + vec1[1]*vec2[1]
					mag1 = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1])
					mag2 = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1])
					#print dotProduct, mag1, mag2, mag1*mag2, dotProduct/(mag1*mag2)
					#cornerAngle 
					try:
						cornerAngle = acos(dotProduct/(mag1*mag2))
	
						len1 = sqrt((p1_1[0]-p1_2[0])**2 + (p1_1[1]-p1_2[1])**2)
						len2 = sqrt((p2_1[0]-p2_2[0])**2 + (p2_1[1]-p2_2[1])**2)
						#if len1 >= 20.0 and len2 >= 20.0:
						#if cornerAngle > 0.25 and cornerAngle < 2.5:
						if True:	
							toDraw.append(line1)
							toDraw.append(line2)
							print i, len1, len2, dist1, dist2, dist3, cornerAngle
							points.append(pnt)
	
					except:
						pass
	
	
		#for line in lines:
		"""
		#for line in toDraw:
		#for line in newLines:
			#print "line:", dist, line[0], line[1]
	
		#for pnt in points:
			#print pnt
			#cv.Circle(color_dst, pnt, 5, cv.CV_RGB(0,0,255), thickness=1, lineType=8, shift=0) 
	
	
		"""
	
		cv.PreCornerDetect(img0, corners, apertureSize=9)
		#cv.SaveImage("%04u_5_2.png" % i, corners)
	
		count = 0
		maxVal = 400
		for m in range(maxX):
			for n in range(maxY):
				if corners[n,m] > maxVal:
					cv.Circle(color_dst, (m,n), 5, cv.CV_RGB(0,0,255), thickness=1, lineType=8, shift=0) 
					count += 1
		print count, "points over", maxVal
	
		filteredPoints = findCornerCandidates(polarLines, edgeCopy, corners)
		#filteredPoints = filterIntersections(polarLines, maxX, maxY)
		#print len(polarLines), "lines produced", len(filteredPoints), "corner points"
	
		for pnt, line1, line2 in filteredPoints:
			indX = cv.Floor(pnt[0])
			indY = cv.Floor(pnt[1])
			#cv.Circle(color_dst, (indX,indY), 5, cv.CV_RGB(0,0,255), thickness=1, lineType=8, shift=0) 
	
			rho1 = line1[0]
			rho2 = line2[0]
			theta1 = line1[1]
			theta2 = line2[1]
			if theta1 > theta2:
				cornerAngle = theta1 - theta2
			else:
				cornerAngle = theta2 - theta1
	
			a = cos(theta1)
			b = sin(theta1)
			x0 = a * rho1 
			y0 = b * rho1
			pt1 = (cv.Round(x0 + 400*(-b)), cv.Round(y0 + 400*(a)))
			pt2 = (cv.Round(x0 - 400*(-b)), cv.Round(y0 - 400*(a)))
			seg1 = (pt1,pt2)
	
			a = cos(theta2)
			b = sin(theta2)
			x0 = a * rho2 
			y0 = b * rho2
			pt1 = (cv.Round(x0 + 400*(-b)), cv.Round(y0 + 400*(a)))
			pt2 = (cv.Round(x0 - 400*(-b)), cv.Round(y0 - 400*(a)))
			seg2 = (pt1,pt2)
	
			#cv.Line(color_dst, seg1[0], seg1[1], cv.CV_RGB(255, 0, 0), 1, 8)
			#cv.Line(color_dst, seg2[0], seg2[1], cv.CV_RGB(255, 0, 0), 1, 8)
			#print "cornerAngle =", cornerAngle
	
		cv.SaveImage("%04u_5.png" % i, color_dst) 
	
		#for line in filteredLines:
		#	cv.Line(color_dst, line[0], line[1], cv.CV_RGB(255, 0, 0), 1, 8)
		#cv.SaveImage("%04u_5_1.png" % i, color_dst) 
	
		#cv.Normalize(corners, corners, 255.0, 0.0)
		#cv.SaveImage("%04u_5_2.png" % i, corners)
	
		# assume that the image is floating-point 
		#corners = cv.CloneMat(image)
		#cv.PreCornerDetect(image, corners, 3)
	
		#dilated_corners = cv.CloneMat(image)
		dilated_corners = cv.CreateImage(cv.GetSize(img0), cv.IPL_DEPTH_32F, 1)
		cv.Dilate(corners, dilated_corners, None, 1)
		#cv.SaveImage("%04u_5_3.png" % i, dilated_corners)
	
		#for m in range(maxX):
		#	for n in range(maxY):
		#		val1, val2 = corners[n,m], dilated_corners[n,m]
		#		if val1 > 0 or val2 > 0:
		#			print n,m,val1, val2
	
		corner_mask = cv.CreateImage(cv.GetSize(img0), cv.IPL_DEPTH_8U, 1)
		#cv.Sub(corners, dilated_corners, corners)
		#cv.CmpS(corners, 0, corner_mask, cv.CV_CMP_GE)
		#cv.SaveImage("%04u_5_4.png" % i, corner_mask)
	
		#corner_mask = cv.CreateMat(size[0], size[1], cv.CV_8UC1)
	
		maxVals = []
		maxLocs = []
		for j in range(1):
			#temp_img = cv.CloneImage(dilated_corners)
			minVal, maxVal, minLoc, maxLoc = cv.MinMaxLoc(dilated_corners)
			#minVal, maxVal, minLoc, maxLoc = cv.MinMaxLoc(temp_img)
			maxVals.append(maxVal)
			maxLocs.append(maxLoc)
	 
			" inhibit location "
			for k in range(-4,4):
				for l in range(-4,4):
					#cv.Set2D(dilated_corners, maxLoc[0] + k, maxLoc[1] + l, 0) 
					#print "setting", (maxLoc[0] + k, maxLoc[1]+l)
					
					#print dilated_corners[maxLoc[0] + k, maxLoc[1] + l] 
					dilated_corners[maxLoc[0] + k, maxLoc[1] + l] = 0
					#print dilated_corners[maxLoc[0] + k, maxLoc[1] + l] 
	
	
		#cv.Circle(color_dst, maxLoc, 5, cv.CV_RGB(0,0,255), thickness=3, lineType=8, shift=0) 
	
		#print maxVals
		#print maxLocs
		#cv.Sub(corners, dilated_corners, corners)
		#cv.CmpS(corners, 0, corner_mask, cv.CV_CMP_GE)
		#corners
		#corner_mask
	
	
		# Find up to 300 corners using Harris
		#for (x,y) in cv.GoodFeaturesToTrack(img0, eig_image, temp_image, 10, 1.0, 1.0, useHarris = True):
		#print cv.GoodFeaturesToTrack(img0, eig_image, temp_image, 10, 1.0, 0.01)
		#for (x,y) in cv.GoodFeaturesToTrack(img0, eig_image, temp_image, 10, 1.0, 1.0):
		#	print i, ":good feature at", x,y
		
	
		#cv.SaveImage("%04u_1.png" % i, img0)
		#cv.SaveImage("%04u_2.png" % i, corners)
		#cv.SaveImage("%04u_2.png" % i, dilated_corners)
		#cv.SaveImage("%04u_3.png" % i, harris_dst)
		#cv.SaveImage("%04u_4.png" % i, edges)
		#cv.SaveImage("%04u_5.png" % i, color_dst) 
		"""
	
	
	#print "channels:", img0.channels
	#print "depth:", img0.depth
	#print "height:", img0.height
	#print "nChannels:", img0.nChannels
	#print "origin:", img0.origin
	#print "width:", img0.width

