cdef extern from "math.h":
	double fabs(double)

cdef inline double lim( double value, double low, double high):
	if value > high:
		return high

	if value < low:
		return low

	return value
	
# input:  self.phi, WORLD_STEP, self.goal_phi, self.init_phi, self.tol, self.pgain, self.igain, self.dgain, self.last_err, self.maxVel
# output: self.errSum, sefl.lastErr, self._done, vel

def commandServo(Phi, step, goalPhi, initPhi, tolerance, PGain, IGain, DGain, lastErr, maxVelocity):	

	cdef double errSum = 0.0
	cdef double vel = 0.0
	cdef double pterm, iterm, dterm, new_last_err
	cdef int isDone = 0
	
	# input arguments
	cdef double pgain = PGain
	cdef double igain = IGain
	cdef double dgain = DGain
	cdef double phi = Phi
	cdef double dt = step
	cdef double goal_phi = goalPhi
	cdef double init_phi = initPhi
	cdef double tol = tolerance
	cdef double last_err = lastErr
	cdef double maxVel = maxVelocity
	
	# compute the servo joint error
	cdef double err = goal_phi - (phi + init_phi)

	if fabs(err) > tol:
		errSum += dt*err
		pterm = pgain*err
		iterm = igain*errSum
		dterm = dgain*(err-last_err)/dt
		new_last_err = err
		vel = lim(pterm+iterm+dterm, -maxVel, maxVel)

	else:
		isDone=True
		errSum=0.0
		new_last_err=0.0
		
	return errSum, new_last_err, isDone, vel

	