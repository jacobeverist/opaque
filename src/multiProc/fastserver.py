# This is the Twisted Fast Poetry Server, version 1.0

import optparse, os

from twisted.internet.protocol import ServerFactory, Protocol


def parse_args():
    usage = """usage: %prog [options] poetry-file

This is the Fast Poetry Server, Twisted edition.
Run it like this:

  python fastpoetry.py <path-to-poetry-file>

If you are in the base directory of the twisted-intro package,
you could run it like this:

  python twisted-server-1/fastpoetry.py poetry/ecstasy.txt

to serve up John Donne's Ecstasy, which I know you want to do.
"""

    parser = optparse.OptionParser(usage)

    help = "The port to listen on. Default to a random available port."
    parser.add_option('--port', type='int', help=help)

    help = "The interface to listen on. Default is localhost."
    parser.add_option('--iface', help=help, default='localhost')

    options, args = parser.parse_args()

    #if len(args) != 1:
    #    parser.error('Provide exactly one poetry file.')

    #poetry_file = args[0]

    #if not os.path.exists(args[0]):
    #    parser.error('No such file: %s' % poetry_file)

    return options


"""
Controller.py
Starts and manages solvers in separate processes for parallel processing.
Provides an interface to the Flex UI.
"""
startIP = 8800
PORT = 8992
import os, sys
from subprocess import Popen
from twisted.spread import pb
from twisted.internet import reactor, defer
from twisted.web import server, resource
#from pyamf.remoting.gateway.twisted import TwistedGateway

class Controller(object):

    # Utilities:
    def broadcastCommand(self, remoteMethodName, arguments, nextStep, failureMessage):
        "Send a command with arguments to all solvers"
        print "broadcasting ...",
        deferreds = [solver.callRemote(remoteMethodName, arguments) for solver in self.solvers.values()]
        print "broadcasted"
        reactor.callLater(3, self.checkStatus)
        # Use a barrier to wait for all to finish before nextStep:
        defer.DeferredList(deferreds, consumeErrors=True).addCallbacks(nextStep,
            self.failed, errbackArgs=(failureMessage))

    def scheduleJobs(self, remoteMethodName, nextStep):
        deferreds = []

	self.jobFinishedCallback = nextStep

	solvers = self.solvers.values()
	
	#for arguments in self.constraintData:
        for i in range(len(self.constraintData)):
            arguments = self.constraintData[i]
            solver = solvers[ i % len(solvers) ] 
            deferreds.append(solver.callRemote(remoteMethodName, arguments))

        defer.DeferredList(deferreds, consumeErrors=True).addCallbacks(self.collectResults,
            self.failed, errbackArgs=("shapeICP call failed"))

    def checkStatus(self):
        "Show that solvers can still receive messages"
        for solver in self.solvers.values():
            solver.callRemote("status").addCallbacks(lambda r: sys.stdout.write(r + "\n"),
                self.failed, errbackArgs=("Status Check Failed"))
        print "Status calls made"

    def failed(self, results, failureMessage="Call Failed"):
        for (success, returnValue), (address, port) in zip(results, self.solvers):
            if not success:
                raise Exception("address: %s port: %d %s" % (address, port, failureMessage))

    def __init__(self):
        self.constraintData = []

        #f = open("shapeInitData.txt")
        #self.constraintData = eval(f.read().rstrip())
        #f.close()
	#self.probIndex = 0
	
        cores = detectCPUs()
        print "Cores:", cores
        # Solver connections will be indexed by (ip, port):
        self.solvers = dict.fromkeys([("localhost", i) for i in range(startIP, startIP + cores)])
        # Start a subprocess on a core for each solver:
        self.pids = [Popen(["python", "solver.py", str(port)]).pid for ip, port in self.solvers]
        print "PIDs:", self.pids
        self.connected = False
        reactor.callLater(1, self.connect) # Give the solvers time to start

    def connect(self):
        "Begin the connection process"
        connections = []
        for address, port in self.solvers:
            factory = pb.PBClientFactory()
            reactor.connectTCP(address, port, factory)
            connections.append(factory.getRootObject())
        defer.DeferredList(connections, consumeErrors=True).addCallbacks(
            self.storeConnections, self.failed, errbackArgs=("Failed to Connect"))

    def storeConnections(self, results):
        for (success, solver), (address, port) in zip(results, self.solvers):
            self.solvers[address, port] = solver
        print "Connected; self.solvers:", self.solvers
        self.connected = True

    def start(self, constraintData, nextStep ):
        self.constraintData = constraintData

        print "starting:", len(self.constraintData)
	#print self.constraintData

        "Begin the solving process"
        if not self.connected:
            return reactor.callLater(0.5, self.start)

        self.scheduleJobs("shapeICP2", nextStep)
        #self.scheduleJobs("shapeICP", nextStep)
        #self.scheduleJobs("shapeICP", self.constraintData[self.probIndex])

    def collectResults(self, results):
        print "received", len(results), "results from", len(self.constraintData), "problems"
        self.jobFinishedCallback(results)
        #reactor.callLater(0, self.jobFinishedCallback) # Give the solvers time to start
        #for i in range(len(results)):
        #     print results[i][0]
        #print "step 3 results:", results

def detectCPUs():
    """
    Detects the number of CPUs on a system. Cribbed from pp.
    """
    # Linux, Unix and MacOS:
    if hasattr(os, "sysconf"):
        if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"):
            # Linux & Unix:
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
        else: # OSX:
            return int(os.popen2("sysctl -n hw.ncpu")[1].read())
    # Windows:
    if os.environ.has_key("NUMBER_OF_PROCESSORS"):
            ncpus = int(os.environ["NUMBER_OF_PROCESSORS"]);
            if ncpus > 0:
                return ncpus
    return 1 # Default

class PoetryProtocol(Protocol):

    def connectionMade(self):
        self.rawData = ''

        #self.transport.write(self.factory.poem)
        #self.transport.loseConnection()

    def connectionLost(self, reason):
        pass
        #print "lost connection:", reason


    def dataReceived(self, data):
        self.rawData += data
        #print "received:", data

	if data[-1] == '\n':
             self.factory.controller.start(eval(self.rawData), self.returnResult)

    def returnResult(self, results):
	
        finalResults = []
        for result in results:
	    if result[0]:
                finalResults.append(result)
            else:
                finalResults.append((False, ([0.0,0.0,0.0], 1e100)))

	self.transport.write(repr(finalResults))
	self.transport.loseConnection()


class PoetryFactory(ServerFactory):

    protocol = PoetryProtocol

    def __init__(self, controller):
        self.poem = ""
        self.controller = controller


def main():
    options = parse_args()

    controller = Controller()

    factory = PoetryFactory(controller)


    from twisted.internet import reactor


    port = reactor.listenTCP(options.port or 10000, factory,
                             interface=options.iface)

    reactor.run()


if __name__ == '__main__':
    main()
