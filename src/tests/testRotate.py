from Image import *
import ImageChops


newImage = open("localOccMap001_0000.png")

for i in range(10):
	degrees = i*10.0

	temp = ImageChops.offset(newImage, int(degrees))

	#print "rotating by", degrees
	#temp = temp.rotate(degrees, BICUBIC)
	temp.save("rotate%02u.png" % i)

