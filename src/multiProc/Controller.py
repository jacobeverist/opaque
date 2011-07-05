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

        f = open("shapeInitData.txt")
        self.constraintData = eval(f.read().rstrip())
        f.close()
	self.probIndex = 0
	
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

    def start(self, results = []):

        "Begin the solving process"
        if not self.connected:
            return reactor.callLater(0.5, self.start)

        self.scheduleJobs("shapeICP", self.constraintData[self.probIndex])

    def collectResults(self, results):
        print "received", len(results), "results from", len(self.constraintData), "problems"
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

class ControlInterface(pb.Root):
    """
    Interface to Flex control panel (Make sure you have at least PyAMF 0.3.1)
    """
    def __init__(self, controller):
        self.controller = controller

    def remote_start(self, _):
        self.controller.start()
        return "Starting parallel jobs"

    def remote_terminate(self, _):
        for solver in controller.solvers.values():
            solver.callRemote("terminate").addErrback(self.controller.failed, "Termination Failed")
        reactor.callLater(1, reactor.stop)
        return "Terminating remote solvers"

if __name__ == "__main__":

    controller = Controller()

    reactor.listenTCP(PORT, pb.PBServerFactory(ControlInterface(controller)))

    # Place the namespace mapping into a TwistedGateway:
    #gateway = TwistedGateway({ "controller": FlexInterface(controller) })
    # Publish the PyAMF gateway at the root URL:
    #root = resource.Resource()
    #root.putChild("", gateway)
    # Tell the twisted reactor to listen:
    #reactor.listenTCP(FlexControlPanelPort, server.Site(root))
    # One reactor runs all servers and clients:
    #controller.start()
    reactor.run()

