# configuration... constants...
import ogre.renderer.OGRE as ogre
import ogre.physics.OgreOde as OgreOde


# python-ogre hacks
#ogre.Vector3().ZERO = ogre.Vector3(0,0,0)
#ogre.Vector3().UNIT_X = ogre.Vector3(1,0,0)
#ogre.Vector3().UNIT_Y = ogre.Vector3(0,1,0)
#ogre.Vector3().UNIT_Z = ogre.Vector3(0,0,1)


# Map Space Parameters
PIXELSIZE = 0.05
#MAPSIZE = 30.0
MAPSIZE = 22.0

# Shadow Mapping parameters
DISTAL_RANGE = 2
INC_VAL = 0.12

# Snake Probe Parameters
ARMLENGTH = 0.15
ARMWIDTH = 0.05
STD_WEIGHT = 1.11
#STD_WEIGHT = 0.1
NUM_SEGS = 40


# Simulation Parameters
#WORLD_STEP = 1e-3
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
SERVO_VMAX = ogre.Math.AngleUnitsToRadians(116.0) # degrees / sec

# servo PID gains
SERVO_PGAIN = 100
SERVO_IGAIN = 1
SERVO_DGAIN = .001

# servo error tolerance
SERVO_ERRMIN = ogre.Math.AngleUnitsToRadians(0.01) # _deg

# servo stops
HI_SERVO_STOP = ogre.Math.AngleUnitsToRadians(160.0) #_deg
LO_SERVO_STOP = ogre.Math.AngleUnitsToRadians(-160.0) #_deg
