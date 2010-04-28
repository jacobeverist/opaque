# configuration... constants...
from math import pi

DEGREES_TO_RADIANS = pi / 180.0


# Map Space Parameters
PIXELSIZE = 0.05
MAPSIZE = 22.0

# Shadow Mapping parameters
DISTAL_RANGE = 2
INC_VAL = 0.12

# Snake Probe Parameters
ARMLENGTH = 0.15
ARMWIDTH = 0.05
STD_WEIGHT = 1.11
NUM_SEGS = 40


# Simulation Parameters
WORLD_STEP = 0.001

ANY_QUERY_MASK                  = 1<<0
ZOMBIE_QUERY_MASK               = 1<<1
GEOMETRY_QUERY_MASK             = 1<<2
VEHICLE_QUERY_MASK              = 1<<3
STATIC_GEOMETRY_QUERY_MASK      = 1<<4

# Servo Parameters

'''
ODE Hinge Joint Parameter Names
    maxForce, maxVelocity, Tol, HiStop, LoStop, PGain, IGain, DGain
'''

# servo force/speed limits
SERVO_FMAX = 2.5 # Newton metres
SERVO_VMAX = 116.0 * DEGREES_TO_RADIANS # degrees / sec

# servo PID gains
SERVO_PGAIN = 100
SERVO_IGAIN = 1
SERVO_DGAIN = .001

# servo error tolerance
SERVO_ERRMIN = 0.01 * DEGREES_TO_RADIANS # _deg

# servo stops
HI_SERVO_STOP = 160.0 * DEGREES_TO_RADIANS # _deg
LO_SERVO_STOP = -160.0 * DEGREES_TO_RADIANS # _deg
