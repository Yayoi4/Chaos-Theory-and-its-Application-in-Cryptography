from vpython import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

canvas(title = "Lorentz Attractor")

#graph for distance graph
fig, ax = plt.subplots()
axs, ays = [],[]

dt = 0.0001 #step size for the calculations
ddt = 0.001 #step size for the rendering (visual)
Pr = 10 #Prandtl number
Ra = 28 #Rayleigh number
Beta = 8/3 #Ratio

Start = False #Boolean that defines if program runs
Exit = False #Boolean that defines if program exits
Started = False #Boolean to check if program already began once (whether the attractor has already been initialized or not)

wtext(text="Number of objects: ") #Text for nc_input
nc_input = winput(bind=None, text="2") #Input box for number of copies, default value 2
wtext(text="\nInitial position: ") #Text for spos_input
spos_input = winput(bind=None, text="1,1,1") #Input box for initial position of the first object, default value (1,1,1)
wtext(text="\nDifference in position (x,y,z) between the copies: ") #Text for d
d_input = winput(bind=None, text="0.01,0,0") #Position distance between the different copies, default value (0.01,0,0)
wtext(text="\n")

class Lorentz:
    def __init__(self, pos, color): #sets up the objects
        self.s = sphere(pos=pos, make_trail=False, trail_type="curve", trail_radius=0.05, max_trail_points=1000000, radius = 0.05, color = color) #creating the object (sphere)
        self.pos  = self.s.pos #copy of positions, that lets the program keep calculating without rendering each step
        self.t = 0 #time counter

    def rk4(self): #Runge-Kutta 4 algorithm
        global Pr,Ra,Beta,dt

        #k1
        #Velocities (x,y,z) for k1
        vx1 = (Pr*(self.pos.y-self.pos.x)) 
        vy1 = (self.pos.x*(Ra-self.pos.z)-self.pos.y)
        vz1 = (self.pos.x*self.pos.y-Beta*self.pos.z)

        #k2
        #Updated positions (x,y,z) at dt/2 using the velocities calculated at k1
        x2 = self.pos.x+vx1*dt/2
        y2 = self.pos.y+vy1*dt/2
        z2 = self.pos.z+vz1*dt/2
        #Velocities (x,y,z) for k2 using the updated positions
        vx2 = (Pr*(y2-x2))
        vy2 = (x2*(Ra-z2)-y2)
        vz2 = (x2*y2-Beta*z2)

        #k3
        #Updated positions (x,y,z) at dt/2 using the velocities calculated at k2
        x3 = self.pos.x+vx2*dt/2
        y3 = self.pos.y+vy2*dt/2
        z3 = self.pos.z+vz2*dt/2
        #Velocities (x,y,z) for k3 using the updated positions
        vx3 = (Pr*(y3-x3))
        vy3 = (x3*(Ra-z3)-y3)
        vz3 = (x3*y3-Beta*z3)

        #k4
        #Updated positions (x,y,z) at dt using the velocities calculated at k3
        x4 = self.pos.x+vx3*dt
        y4 = self.pos.y+vy3*dt
        z4 = self.pos.z+vz3*dt
        #Velocities (x,y,z) for k4 using the updated positions
        vx4 = (Pr*(y4-x4))
        vy4 = (x4*(Ra-z4)-y4)
        vz4 = (x4*y4-Beta*z4)

        #calculating the final velocity (x,y,z)
        vx = (vx1+2*vx2+2*vx3+vx4)/6
        vy = (vy1+2*vy2+2*vy3+vy4)/6
        vz = (vz1+2*vz2+2*vz3+vz4)/6

        #calculate the final updated positions with that final velocity
        self.pos.x += vx*dt
        self.pos.y += vy*dt
        self.pos.z += vz*dt

        #update the counter
        self.t += dt
    
    def update(self):
        for i in range(int(1/dt*ddt)): #calculating steps between every visual update
            self.rk4()
        self.s.pos = self.pos #updating the visual
        if (self.t > 500*dt): #late start of the trail for smoothness. 
            self.s.make_trail=True

bodies = [] #list of all copies
colors = [color.blue, color.red, color.green, color.yellow] #list of colors for every copy
def StartSim():
    global Start, Exit, nc, bodies, colors, Started
    if not Started:
        spos_list = list(map(float,spos_input.text.split(","))) #change text from the initial position input box to float numbers
        spos = vector(spos_list[0],spos_list[1],spos_list[2]) #define the starting position vector using spos_list
        d_list = list(map(float,d_input.text.split(","))) #change text from d_input to float numbers
        d = vector(d_list[0],d_list[1],d_list[2]) #define difference vector from d_list
        for i in range(int(nc_input.text)): #create all objects using the inputs given
            bodies.append(Lorentz(spos + i*d,colors[i]))
        Started = True #Attractor has been initialized
    print("Start Simulation")
    Start = True #begin simulation

def ExitProgram(): #Exit program
    global Exit
    Exit = True
    Start = False
    print("Ending Simulation")

def pause(): #Pause the program
    global Start
    Start = False
    print("Pause...")

StartButton = button(bind=StartSim, text='Start Simulation!') #Button to start the simulation
EndButton = button(bind=ExitProgram, text="Exit!") #Button to exit the simulation
Pause = button(bind=pause, text="Pause") #Button to pause the simulation

def dist(a,b):
    return sqrt(pow(abs(a.s.pos.x-b.s.pos.x),2)+pow(abs(a.s.pos.y-b.s.pos.y),2)+pow(abs(a.s.pos.z-b.s.pos.z),2))

while True:
    rate(1/ddt/5) #How often it loops per second, limited because of the limited speed at which the computer can render the visual
    if Start == True: #when Start button is pressed
        for b in bodies:
            b.update()#update every object
        if int(nc_input.text)>=2: #if at least 2 copies, plot graph showing the position difference between the s1 of the first two copies
            axs.append(bodies[0].t)
            ays.append(dist(bodies[0],bodies[1]))
    if Exit == True: #Exit the Program
        break

ax.plot(axs,ays)
ax.set(xlabel="t", ylabel="Distance")
fig.savefig("Distance_LA.pdf") #save plot as pdf
plt.show() #show the matplotlib plots