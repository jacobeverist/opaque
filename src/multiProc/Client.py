"Client.py: Uses the calculation service across the network"
import sys
from twisted.spread import pb
from twisted.internet import reactor, defer

PORT = 8992
COUNT = 3

class ModelCalculator:

	def __init__(self, host):
		factory = pb.PBClientFactory()
		reactor.connectTCP(host, PORT, factory)
		factory.getRootObject().addCallbacks(self.connected, self.failure)

	def connected(self, perspective):
		perspective.callRemote('start', 0).addCallbacks(self.success, self.failure)
		#self.deferred = perspective.callRemote('start', 0).addCallbacks(self.success, self.failure)
		#self.deferred.addCallbacks(self.success,self.failure)

		print "connected"

	def success(self, result):
		print "success"
		print result
		reactor.stop()

	def failure(self, _):
		print "remote failure"
		reactor.stop()

ModelCalculator("127.0.0.1")
reactor.run()

#reload(pb)
#reload(reactor)
#reload(defer)

#ModelCalculator("127.0.0.1")
#reactor.run()

print sys.getrefcount(pb)
print sys.getrefcount(reactor)
print sys.getrefcount(defer)

del ModelCalculator
del sys.modules["twisted"]
del sys.modules["twisted.internet"]
del sys.modules["twisted.spread"]

#del sys.modules["pb"]
#del sys.modules["reactor"]
#del sys.modules["defer"]

print sys.getrefcount(pb)
print sys.getrefcount(reactor)
print sys.getrefcount(defer)


