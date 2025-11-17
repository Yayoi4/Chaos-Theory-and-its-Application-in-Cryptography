from vpython import *
import numpy as np

canvas(title="Billiards")

# Parameters
dx = 4
dy = 4
r = 0.05
vx = 2
vy = 3
sx = 0
sy = 0
dt = 0.001
t = 0
Start = False
Exit = False

class billiards:
    def __init__(self, pos=vector(sx,sy,0), vel=vector(vx,vy,0), colors = color.green):
        self.ball = sphere(pos=pos, vel=vel, radius=r, color=colors)
        attach_trail(self.ball, retain=100000, color=colors)
        self.pos = pos
        self.vel = vel

    def CheckWalls(self): #check if there's a collision against the wall
        pos = self.ball.pos
        vel = self.ball.vel

        # wall collision - Y
        if abs(pos.y) >= dy/2 - r:
            vel.y = -vel.y
            #correct path if over the border
            if pos.y > dy/2 - r:
                pos.y = dy/2 - r - (pos.y - (dy/2 - r))
            elif pos.y < -dy/2 + r:
                pos.y = -dy/2 + r - (pos.y + dy/2 - r)

        # wall collision - X
        if abs(pos.x) >= dx/2 - r:
            vel.x = -vel.x
            #correct path if over the border
            if pos.x > dx/2 - r:
                pos.x = dx/2 - r - (pos.x - (dx/2 - r))
            elif pos.x < -dx/2 + r:
                pos.x = -dx/2 + r - (pos.x + dx/2 - r)

        self.ball.pos = pos
        self.ball.vel = vel

    def update(self):
        self.ball.pos = self.ball.pos + self.ball.vel * dt
        self.CheckWalls()
        

balls = []

def StartSimulation():
    global Start, balls
    lw = box(pos=vector(-dx/2,0,0), length=0.1, height=dy+0.1, width=0.1, color=color.blue)
    tw = box(pos=vector(0,dy/2,0), length=dx+0.1, height=0.1, width=0.1, color=color.blue)
    rw = box(pos=vector(dx/2,0,0), length=0.1, height=dy+0.1, width=0.1, color=color.blue)
    bw = box(pos=vector(0,-dy/2,0), length=dx+0.1, height=0.1, width=0.1, color=color.blue)
    balls = [
        billiards(vector(0.02, 0.0, 0), vector(2, 1.5, 0)),
        billiards(vector(0.01, 0, 0), vector(2, 1.5, 0)),
        billiards(vector(0, 0, 0), vector(2, 1.5, 0))
    ]
    Start = True

def ExitSimulation():
    global Exit
    Exit = True

def PauseSimulation():
    global Start
    Start = False

Start = button(bind=StartSimulation, text="Start")
Pause = button(bind=PauseSimulation, text="Pause")
Exist = button(bind=ExitSimulation, text="Exit")

# Simulation
while True:
    rate(300)
    # pairwise collisions
    if Start == True:
        for b in balls:
            b.update()
        t += dt
    if Exit == True:
        break


