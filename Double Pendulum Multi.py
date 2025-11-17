from vpython import *
import numpy as np
import sys
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib as mpl

canvas(title='Double Pendulum')

fig, ax = plt.subplots() #Total Energy Plot
#fig, bx = plt.subplots()
#fig, cx = plt.subplots()
fig, dx = plt.subplots()
dxs, dys = [],[]
axs, ays = [],[]
bxs, bys = [],[]
cxs, cys = [],[]


g = 9.81 #gravity constant
dt = 0.0001 #step size for the calculations
ddt = 0.01 #step size for the rendering (visual)
t = 0 #time counter

Start = False #Boolean that defines if program runs
Exit = False #Boolean that defines if program exits
Started = False #Boolean to check if program already began once (whether the pendulum has already been initialized or not)

wtext(text="\nNumber of copies: ") #Text for nc_input
nc_input = winput(bind=None, text="2") #Input box for number of copies, default value 2
wtext(text="\nMasses of the two balls: ") #Text for m1_input&m2_input
m1_input = winput(bind=None, text="1.0") #Mass of the first ball
m2_input = winput(bind=None, text="1.0") #Mass of the second ball
wtext(text="\nStarting angle of the balls: ") #Text for theta1_input&theta2_input
theta1_input = winput(bind=None, text="90") #Starting angle for the first ball
theta2_input = winput(bind=None, text="90") #Starting angle for the second ball
wtext(text="\nLength of the rods: ") #Text for l1_input&l2_input
l1_input = winput(bind=None, text="1.0") #Length of the first rod
l2_input = winput(bind=None, text="1.0") #Length of the second rod
wtext(text="\nDifference in masses for the two balls, lengths of the two rods, and starting angles of the balls between the copies: ") #Text for d_input
d_input = winput(bind=None, text="0,0,0,0,0,0.01") #Difference between the copies for the masses, lengths and starting angles
wtext(text="\n")
gr1 = gcurve(visible = True) #Graph for Total Energy
#gr2 = gcurve(visible = True)
wtext(text="\n")

class DoublePendulum:
    def __init__(self, m1=1, m2=1, L1=1, L2=1, theta1=90*pi/180, theta2=90*pi/180, omega1=0, omega2=0, origin=vector(0,0,0), color=color.red):
        #initiate all the initial conditions
        self.m1 = m1
        self.m2 = m2
        self.L1 = L1
        self.L2 = L2
        self.theta1 = theta1
        self.theta2 = theta2
        self.omega1 = omega1
        self.omega2 = omega2
        self.origin = origin

        #create the pendulum (root, balls, rods)
        self.s1 = sphere(radius=0.1, color=color)
        self.s2 = sphere(radius=0.1, color=color, make_trail=False, max_trail_points = 1000)
        self.root = sphere(pos=origin, radius=0.1, color=color)
        self.pendulum1 = cylinder(radius=0.015, color=color)
        self.pendulum2 = cylinder(radius=0.015, color=color)

        #initialize the initial position
        self.s1.pos = self.root.pos + vector(L1*sin(theta1), -L1*cos(theta1), 0)
        self.s2.pos = self.s1.pos+vector(L2*sin(theta2), -L2*cos(theta2), 0)
        self.pendulum1.axis = self.s1.pos-self.root.pos
        self.pendulum2.pos = self.s1.pos
        self.pendulum2.axis = self.s2.pos-self.s1.pos 

    def f1(self, theta1, theta2, omega1, omega2): #f1 function (mentioned in the paper)
        global g
        return ((-self.L2/self.L1)*(self.m2/(self.m1+self.m2))*pow(omega2,2)*sin(theta1-theta2)-g/self.L1*sin(theta1))

    def f2(self, theta1, theta2, omega1, omega2): #f2 function (mentioned in the paper)
        return ((self.L1/self.L2)*pow(omega1,2)*sin(theta1-theta2)-(g/self.L2)*sin(theta2))

    def alpha1(self, theta1, theta2): #alpha1 function (mentioned in the paper)
        return ((self.L2/self.L1)*(self.m2/(self.m1+self.m2))*cos(theta1-theta2))
    
    def alpha2(self, theta1, theta2): #alpha2 function (mentioned in the paper)
        return ((self.L1/self.L2)*cos(theta1-theta2))
    
    def g1(self, theta1, theta2, omega1, omega2): #g1 function (mentioned in the paper)
        return ((self.f1(theta1, theta2, omega1, omega2)-self.alpha1(theta1, theta2)*self.f2(theta1, theta2, omega1, omega2))/(1-self.alpha1(theta1, theta2)*self.alpha2(theta1, theta2)))

    def g2(self, theta1, theta2, omega1, omega2): #g2 function (mentioned in the paper)
        return ((-self.alpha2(theta1, theta2)*self.f1(theta1, theta2, omega1, omega2)+self.f2(theta1,theta2, omega1, omega2))/(1-self.alpha1(theta1, theta2)*self.alpha2(theta1, theta2)))

    def rk4(self): #Runge-Kutta 4 algorithm
        global dt

        #defining shortcuts
        m1, m2, L1, L2 = self.m1, self.m2, self.L1, self.L2
        theta1, theta2 = self.theta1, self.theta2
        omega1, omega2 = self.omega1, self.omega2

        #rk1
        #define velocity for k1
        aomega1 = omega1
        aomega2 = omega2
        #calculate acceleration for k1
        aomegadot1 = self.g1(theta1, theta2, omega1, omega2)
        aomegadot2 = self.g2(theta1, theta2, omega1, omega2)

        #rk2
        #updated angles and velocities from k1 values
        btheta1 = theta1 + aomega1*dt/2
        btheta2 = theta2 + aomega2*dt/2
        bomega1 = omega1 + aomegadot1*dt/2
        bomega2 = omega2 + aomegadot2*dt/2
        #calcualte acceleration for k2
        bomegadot1 = self.g1(btheta1,btheta2,bomega1, bomega2)
        bomegadot2 = self.g2(btheta1,btheta2,bomega1, bomega2)

        #rk3
        #updated angles and velocities from k2 values
        ctheta1 = theta1 + bomega1*dt/2
        ctheta2 = theta2 + bomega2*dt/2
        comega1 = omega1 + bomegadot1*dt/2
        comega2 = omega2 + bomegadot2*dt/2
        #calculate acceleration for k3
        comegadot1 = self.g1(ctheta1,ctheta2,comega1, comega2)
        comegadot2 = self.g2(ctheta1,ctheta2,comega1, comega2)

        #rk4
        #updated angles and velocities from k3 values
        dtheta1 = theta1 + comega1*dt
        dtheta2 = theta2 + comega2*dt
        domega1 = omega1 + comegadot1*dt
        domega2 = omega2 + comegadot2*dt
        #calculate acceleration for k4
        domegadot1 = self.g1(dtheta1,dtheta2,domega1, domega2)
        domegadot2 = self.g2(dtheta1,dtheta2,domega1, domega2)

        #calculate and update final velocity and acceleration
        self.theta1 += (aomega1+2*bomega1+2*comega1+domega1)/6*dt
        self.theta2 += (aomega2+2*bomega2+2*comega2+domega2)/6*dt
        self.omega1 += (aomegadot1+2*bomegadot1+2*comegadot1+domegadot1)/6*dt
        self.omega2 += (aomegadot2+2*bomegadot2+2*comegadot2+domegadot2)/6*dt

    def euler(self): #Euler algorithm (for comparison)
        #defining shortcuts
        m1, m2, L1, L2 = self.m1, self.m2, self.L1, self.L2
        theta1, theta2 = self.theta1, self.theta2
        omega1, omega2 = self.omega1, self.omega2

        #calculating new acceleration using Euler algorithm with equations of motion (long version)
        omega1dot = (-g*(2*m1+m2)*sin(theta1)-m2*g*sin(theta1-2*theta2)-2*sin(theta1-theta2)*m2*(omega2*omega2*L2+omega1*omega1*L1*cos(theta1-theta2)))/(L1*(2*m1+m2-m2*cos(2*theta1-2*theta2)))
        omega2dot = (2*sin(theta1-theta2)*(omega1*omega1*L1*(m1+m2)+g*(m1+m2)*cos(theta1)+omega2*omega2*L2*m2*cos(theta1-theta2)))/(L2*(2*m1+m2-m2*cos(2*theta1-2*theta2)))

        #update values
        omega1 += omega1dot * dt
        omega2 += omega2dot * dt
        theta1 += omega1 * dt
        theta2 += omega2 * dt

        self.omega1 = omega1
        self.omega2 = omega2
        self.theta1 = theta1
        self.theta2 = theta2

    def update(self):
        global gr1, gr2, ax, axs, ays, bxs, bys, cxs, cys
        
        
        for i in range(int(1/dt*ddt)): #calculating steps between every visual update
            self.rk4()

        #defining shortcuts
        m1, m2, L1, L2 = self.m1, self.m2, self.L1, self.L2
        theta1, theta2 = self.theta1, self.theta2
        omega1, omega2 = self.omega1, self.omega2

        #updating visual (positions)
        self.s1.pos = self.origin + vector(L1 * sin(theta1), -L1 * cos(theta1), 0)
        self.s2.pos = self.s1.pos + vector(L2 * sin(theta2), -L2 * cos(theta2), 0)
        self.pendulum1.axis = self.s1.pos - self.origin
        self.pendulum2.pos = self.s1.pos
        self.pendulum2.axis = self.s2.pos - self.s1.pos

        #calculating position (from bottom) and velocity
        y1 = L1+L2+(-L1 * cos(theta1))
        y2 = L1+L2+(-L1 * cos(theta1) - L2 * cos(theta2))
        v1x = L1 * omega1 * cos(theta1)
        v1y = L1 * omega1 * sin(theta1)
        v2x = L1 * omega1 * cos(theta1) + L2 * omega2 * cos(theta2)
        v2y = L1 * omega1 * sin(theta1) + L2 * omega2 * sin(theta2)
        #calculating kinetic and potential energy
        T = 1/2 * m1 * (pow(v1x,2) + pow(v1y,2)) + 1/2 * m2 * (pow(v2x,2) + pow(v2y,2))
        V = m1*g*y1 + m2*g*y2
        gr1.plot(t,T+V) #Total energy plot (VPython plot)
        axs.append(t) #add time spot for x-axis matplot
        ays.append(T+V) #add total energy for y-axis for matplot
        #gr2.plot(theta1,theta2)
        #bxs.append(theta1)
        #bys.append(v1x)
        #cxs.append(theta2)
        #cys.append(v2x)

        if t>=0.5: #late start of the trail for smoothness
            self.s2.make_trail=True


bodies = [] #list of all copies
colors = [color.blue, color.red, color.green, color.yellow] #list of colors for every copy
def StartSim(): #Start the Simulation
    global Start, Started, Exit, bodies, colors, t
    if not Started:
        d = list(map(float,d_input.text.split(","))) #list of the difference between every copy
        for i in range(int(nc_input.text)): #create all objects using the inputs given
            bodies.append(DoublePendulum(float(m1_input.text)+i*d[0],float(m2_input.text)+i*d[1],float(l1_input.text)+i*d[2],float(l2_input.text)+i*d[3],float(theta1_input.text)*pi/180+i*d[4],float(theta2_input.text)*pi/180+i*d[5],color=colors[i]))
        # if (int(nc_input.text) >= 2):
        #     gr2.visible = True
        Started = True #Pendulum has been initialized
    print("Start Simulation")
    Start = True #Begin simulation

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

wtext(text="\n\n")

def dist(a,b):
    return sqrt(pow(abs(a.s2.pos.x-b.s2.pos.x),2)+pow(abs(a.s2.pos.y-b.s2.pos.y),2))

while True:
    rate(1/ddt)  #How often it loops per second, limited because of the limited speed at which the computer can render the visual
    if Start == True: #when Start button is pressed
        for b in bodies:
            b.update() #update every copy
        t += ddt #update time
        if int(nc_input.text)>=2:
            dxs.append(t)
            dys.append(dist(bodies[0],bodies[1]))
    if Exit: #Exit the program
        break


#bx.plot(bxs,bys)
#bx.margins(y=0.1)
#cx.plot(cxs,cys)
ax.plot(axs,ays) #plot total energy
ax.margins(y=0.1) #add small margin at the borders
ax.set(xlabel="t", ylabel="Total Energy")
dx.plot(dxs,dys)
dx.set(xlabel="t", ylabel="Distance")
fig.savefig("Distance.pdf") #save plot as pdf
plt.show() #show the matplotlib plots

    
