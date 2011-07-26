"""
solver.py
Solves one portion of a problem, in a separate process on a separate CPU
"""
import sys, random, math
from twisted.spread import pb
from twisted.internet import reactor

#import gen_icp
from opaque.maps import gen_icp


class Solver(pb.Root):

    def __init__(self, id):
        self.id = id

    def __str__(self): # String representation
        return "Solver %s" % self.id

    def remote_initialize(self, initArg):
        return "%s initialized" % self

    def remote_shapeICP2(self, args):
        estPose1 = args[0]
	estPose2 = args[1]
	hull1 = args[2]
	hull2 = args[3]
	posture1_unstable = args[4]
	posture1_stable = args[5]
	posture2_unstable = args[6]
	posture2_stable = args[7]
	costThresh = args[8]
	minMatchDist = args[9]
	plotIteration = args[10]
	n1 = args[11]
	n2 = args[12]

	offset, cost, status = gen_icp.shapeICP2(estPose1, estPose2, hull1, hull2, posture1_unstable, posture1_stable, posture2_unstable, posture2_stable, costThresh, minMatchDist, plotIteration, n1, n2)

        print "sensor constraint: %d -> %d" %(n1,n2), "cost =", cost

	" cast the values to float so they pass the type security of jelly "
	offset = [float(offset[0]), float(offset[1]), float(offset[2])]

	return offset, cost, status

    def remote_shapeICP(self, args):
        estPose1 = args[0]
	firstGuess = args[1]
	hull1 = args[2]
	hull2 = args[3]
	stablePoints1 = args[4]
	stablePoints2 = args[5]
	uniform1 = args[6]
	uniform2 = args[7]
	isForward = args[8]
	(headPoint1,tailPoint1,headPoint2,tailPoint2) = args[9]
	circles = args[10]
	costThresh = args[11]
	minMatchDist = args[12]
	plotIteration = args[13]
	n1 = args[14]
	n2 = args[15]
        offset, cost = gen_icp.shapeICP(estPose1, firstGuess, hull1, hull2, stablePoints1, stablePoints2, uniform1, uniform2, isForward, (headPoint1,tailPoint1,headPoint2,tailPoint2), circles, costThresh, minMatchDist, plotIteration, n1, n2)
        print "sensor constraint: %d -> %d" %(n1,n2), "cost =", cost

	" cast the values to float so they pass the type security of jelly "
	offset = [float(offset[0]), float(offset[1]), float(offset[2])]

	return offset, cost

    def step(self, arg):
        "Simulate work and return result"
        result = 0
        for i in range(random.randint(1000000, 3000000)):
            angle = math.radians(random.randint(0, 45))
            result += math.tanh(angle)/math.cosh(angle)
        return "%s, %s, result: %.2f" % (self, str(arg), result)

    # Alias methods, for demonstration version:
    remote_step1 = step
    remote_step2 = step
    remote_step3 = step

    def remote_status(self):
        return "%s operational" % self

    def remote_terminate(self):
        reactor.callLater(0.5, reactor.stop)
        return "%s terminating..." % self

if __name__ == "__main__":
    port = int(sys.argv[1])
    reactor.listenTCP(port, pb.PBServerFactory(Solver(sys.argv[1])))
    reactor.run()

