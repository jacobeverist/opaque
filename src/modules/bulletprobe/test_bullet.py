
from bulletsnakeprobe import BulletSnakeProbe


quat_param = [0.0, 1.0, 0.0, 0.0]
pos = [1.65, 1.04, 0.0]


bullet_snake = BulletSnakeProbe(quat_param, pos, 40, 0.15, 0.2, 0.15, 0.7)

for i in range(300):
	bullet_snake.frameStarted()

	val = bullet_snake.getGlobalPosition(9)
	print val

exit()

val = bullet_snake.getGlobalOrientation(9)
print val

print bullet_snake.getServo(9)

bullet_snake.frameStarted()

bullet_snake.setServo(9,0.0)

print bullet_snake.getServo(9)

bullet_snake.frameStarted()

print bullet_snake.getServo(9)


