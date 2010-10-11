
from nxsnakeprobe import NxSnakeProbe

nx_snake = NxSnakeProbe()
nx_snake.frameStarted()

val = nx_snake.getGlobalPosition(9)
print val

val = nx_snake.getGlobalOrientation(9)
print val

print nx_snake.getServo(9)

nx_snake.frameStarted()

nx_snake.setServo(9,0.0)

print nx_snake.getServo(9)

nx_snake.frameStarted()

print nx_snake.getServo(9)


