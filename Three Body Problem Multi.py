from vpython import *
import numpy as np
import sys
import random
import matplotlib.pyplot as plt

canvas(title='Three Body Problem')

#graph for distance graph
fig, ax = plt.subplots()
axs, ays = [],[]

G = 1.0 #gravitational constant
dt = 0.0001 #step size for the calculations
ddt = 0.01 #step size for the rendering (visual)

Start = False #Boolean that defines if program runs
Exit = False #Boolean that defines if program exits
Started = False #Boolean to check if program already began once (whether the objects have already been initialized or not)

#Stable solutions
IA1 = [[vector(-1.0024277970,  0.0041695061, 0.0), vector(1.0024277970, -0.0041695061, 0.0), vector(0,0,0)],[vector(0.3489048974 , 0.5306305100, 0.0), vector(0.3489048974, 0.5306305100, 0.0), vector(-2*0.3489048974, -2*0.5306305100, 0.0)]]
IA2 = [[vector(-1.0005576155,  -0.0029240248, 0.0), vector(1.0005576155, 0.0029240248, 0.0), vector(0,0,0)],[vector(0.3064392516 , 0.1263673939, 0.0), vector(0.3064392516, 0.1263673939, 0.0), vector(-2*0.3064392516, -2*0.1263673939, 0.0)]]
IA17 = [[vector(-1.0074958476, 0.0081648176, 0.0), vector(1.0074958476, -0.0081648176, 0.0), vector(0,0,0)],[vector(0.1883232887, 0.5834831526, 0.0), vector(0.1883232887, 0.5834831526, 0.0), vector(-2*0.1883232887, -2*0.5834831526, 0.0)]]
IA77 = [[vector(-1,0,0), vector(1,0,0), vector(0,0,0)],[vector(0.4159559963,0.2988672319,0.0),vector(0.4159559963,0.2988672319,0.0),vector(-2*0.4159559963,-2*0.2988672319,0.0)]]
IA115 = [[vector(-1,0,0), vector(1,0,0), vector(0,0,0)],[vector(0.3369172422,0.2901238678,0.0),vector(0.3369172422,0.2901238678,0.0),vector(-2*0.3369172422,-2*0.2901238678,0.0)]]
IB1 = [[vector(-0.9989071137,  -0.0001484864, 0.0), vector(0.9989071137, 0.0001484864, 0.0), vector(0,0,0)],[vector(0.4646402601  , 0.3963456869, 0.0), vector(0.4646402601, 0.3963456869, 0.0), vector(-2*0.4646402601, -2*0.3963456869, 0.0)]]
IB2 = [[vector(-1.1770534081, -0.5225957568,0.0), vector(1.1770534081, 0.5225957568,0.0), vector(0,0,0)],[vector(0.2446132140,0.3305126876,0.0),vector(0.2446132140,0.3305126876,0.0),vector(-2*0.2446132140,-2*0.3305126876,0.0)]]
IB6 = [[vector(-1.0043366457,  0.0085104316, 0.0), vector(1.0043366457, -0.0085104316, 0.0), vector(0,0,0)],[vector(0.3857847594 , 0.3732858410, 0.0), vector(0.3857847594, 0.3732858410, 0.0), vector(-2*0.3857847594, -2*0.3732858410, 0.0)]]
IIC1 = [[vector(-0.9826146484,  -0.0411837391, 0.0), vector(0.9826146484, 0.0411837391, 0.0), vector(0,0,0)],[vector(0.2710001824 , 0.3415940623, 0.0), vector(0.2710001824, 0.3415940623, 0.0), vector(-2*0.2710001824, -2*0.3415940623, 0.0)]]

#UI
wtext(text="\nNumber of copies: ") #Text for nc_input
nc_input = winput(bind=None, text="2") #Input box for number of copies, dfault value 2
wtext(text="\nInitial positions of the three bodies: ") #Text for pi_input
p1_input = winput(bind=None, text="-1,0,0") #Input box for position of first body
p2_input = winput(bind=None, text="0,1,0") #Input box for position of second body
p3_input = winput(bind=None, text="1,0,0") #Input box for position of third body
wtext(text="\nMasses of the three bodies: ") #Text for mi_input
m1_input = winput(bind=None, text="1.0") #Input box for mass of first body
m2_input = winput(bind=None, text="1.0") #Input box for mass of second body
m3_input = winput(bind=None, text="1.0") #Input box for mass of third body
wtext(text="\nVelocities of the three bodies: ") #Text for vi_input
v1_input = winput(bind=None, text="-1,0,0") #Input box for velocity of first body
v2_input = winput(bind=None, text="0,0,0") #Input box for velocity of second body
v3_input = winput(bind=None, text="1,0,0") #Input box for velocity of third body
wtext(text="\nDifferent in masses of the copies: ")#text for dm_input
dm_input = winput(bind=None, text="0.1,0,0") #Input box for difference in masses of the three bodies between the copies
wtext(text="\nDifferent in position between the copies: ") #text for dpi_input
#Input boxes for difference in position for the three bodies
dp1_input = winput(bind=None, text="0,0,0")
dp2_input = winput(bind=None, text="0,0,0")
dp3_input = winput(bind=None, text="0,0,0")
wtext(text="\nDifferent in velocity between the copies: ") #text for dvi_input
#Input boxes for difference in velocity for the three bodies
dv1_input = winput(bind=None, text="0,0,0")
dv2_input = winput(bind=None, text="0,0,0")
dv3_input = winput(bind=None, text="0,0,0")
wtext(text="\n")

#position and velocity vectors
p1 = vector(0,0,0)
p2 = vector(0,0,0)
p3 = vector(0,0,0)
v1 = vector(0,0,0)
v2 = vector(0,0,0)
v3 = vector(0,0,0)

class ThreeBody:
    def __init__(self, p1, p2, p3, m1=1, m2=1, m3=1, v1=1, v2=1, v3=1, scolor=color.blue):
        #initiate all the initial conditions
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.a1 = vector(0,0,0)
        self.a2 = vector(0,0,0)
        self.a3 = vector(0,0,0)

        self.t = 0 #counter

        #create the three bodies
        self.s1 = sphere(pos=p1, radius = 0.05, color=scolor, s1s2 = vector(0,0,0), s1s3 = vector(0,0,0))
        self.s2 = sphere(pos=p2, radius = 0.05, color=scolor, s2s1 = vector(0,0,0), s2s3 = vector(0,0,0))
        self.s3 = sphere(pos=p3, radius = 0.05, color=scolor, s3s1 = vector(0,0,0), s3s2 = vector(0,0,0))
        attach_trail(self.s1, retain=100000, color=self.s1.color)
        attach_trail(self.s2, retain=100000, color=self.s2.color)
        attach_trail(self.s3, retain=100000, color=self.s3.color)
        #add force vectors between the bodies 
        attach_arrow(self.s1, "s1s2", color=self.s1.color, scale = 0.2, shaftwidth=self.s1.radius/3)
        attach_arrow(self.s1, "s1s3", color=self.s1.color, scale = 0.2, shaftwidth=self.s1.radius/3)
        attach_arrow(self.s2, "s2s1", color=self.s2.color, scale = 0.2, shaftwidth=self.s1.radius/3)
        attach_arrow(self.s2, "s2s3", color=self.s2.color, scale = 0.2, shaftwidth=self.s1.radius/3)
        attach_arrow(self.s3, "s3s1", color=self.s3.color, scale = 0.2, shaftwidth=self.s1.radius/3)
        attach_arrow(self.s3, "s3s2", color=self.s3.color, scale = 0.2, shaftwidth=self.s1.radius/3)

    def dist(self,a,b): #calculate the distance between two bodies
        d1=a.x-b.x
        d2=a.y-b.y
        d3=a.z-b.z
        return sqrt(pow(d1,2) + pow(d2,2) + pow(d3,2))

    def calca(self,ts1,ts2,ts3): #calculate the accelerations of the bodies
        global G
        ta1 = -G*self.m2*(ts1-ts2)/(pow(self.dist(ts1,ts2),3))-G*self.m3*(ts1-ts3)/(pow(self.dist(ts1,ts3),3))
        ta2 = -G*self.m3*(ts2-ts3)/(pow(self.dist(ts2,ts3),3))-G*self.m1*(ts2-ts1)/(pow(self.dist(ts2,ts1),3))
        ta3 = -G*self.m1*(ts3-ts1)/(pow(self.dist(ts3,ts1),3))-G*self.m2*(ts3-ts2)/(pow(self.dist(ts3,ts2),3))
        return ta1, ta2, ta3

    def rk4(self): #Runge-Kutta 4 algorithm
        global G, dt

        #defining shortcuts
        tp1 = self.p1
        tp2 = self.p2
        tp3 = self.p3
        tv1 = self.v1
        tv2 = self.v2
        tv3 = self.v3

        #k1
        #define position and velocity for k1
        as1 = tp1
        as2 = tp2
        as3 = tp3
        av1 = tv1
        av2 = tv2
        av3 = tv3
        #calculate acceleration for k1
        aa1,aa2,aa3 = self.calca(as1,as2,as3)

        #k2
        #update position and veloctiy for k2
        bs1 = tp1+av1*dt/2
        bs2 = tp2+av2*dt/2
        bs3 = tp3+av3*dt/2
        bv1 = tv1+dt/2*aa1
        bv2 = tv2+dt/2*aa2
        bv3 = tv3+dt/2*aa3
        #calculate acceleration for k2
        ba1,ba2,ba3 = self.calca(bs1,bs2,bs3)
        
        #k3
        #update position and velocity for k3
        cs1 = tp1+bv1*dt/2
        cs2 = tp2+bv2*dt/2
        cs3 = tp3+bv3*dt/2
        cv1 = tv1+dt/2*ba1
        cv2 = tv2+dt/2*ba2
        cv3 = tv3+dt/2*ba3
        #calculate acceleration for k3
        ca1,ca2,ca3 = self.calca(cs1,cs2,cs3)

        #k4
        #update position and velocity for k4
        ds1 = tp1+cv1*dt/2
        ds2 = tp2+cv2*dt/2
        ds3 = tp3+cv3*dt/2
        dv1 = tv1+dt/2*ca1
        dv2 = tv2+dt/2*ca2
        dv3 = tv3+dt/2*ca3
        #calculate acceleration for k4
        da1,da2,da3 = self.calca(ds1,ds2,ds3)

        #final velocity
        self.v1 = tv1 + (aa1+2*ba1+2*ca1+da1)/6*dt
        self.v2 = tv2 + (aa2+2*ba2+2*ca2+da2)/6*dt
        self.v3 = tv3 + (aa3+2*ba3+2*ca3+da3)/6*dt

        #final position
        self.p1 = tp1 + (av1+2*bv1+2*cv1+dv1)/6*dt
        self.p2 = tp2 + (av2+2*bv2+2*cv2+dv2)/6*dt
        self.p3 = tp3 + (av3+2*bv3+2*cv3+dv3)/6*dt

    def forceArrow(self): #arrows between bodies to display forces
        global G
        self.s1.s1s2 = -G*self.m2*(self.s1.pos-self.s2.pos)/(pow(self.dist(self.s1.pos,self.s2.pos),3))
        self.s1.s1s3 = -G*self.m3*(self.s1.pos-self.s3.pos)/(pow(self.dist(self.s1.pos,self.s3.pos),3))
        self.s2.s2s1 = -G*self.m1*(self.s2.pos-self.s1.pos)/(pow(self.dist(self.s2.pos,self.s1.pos),3))
        self.s2.s2s3 = -G*self.m3*(self.s2.pos-self.s3.pos)/(pow(self.dist(self.s2.pos,self.s3.pos),3))
        self.s3.s3s1 = -G*self.m3*(self.s3.pos-self.s1.pos)/(pow(self.dist(self.s3.pos,self.s1.pos),3))
        self.s3.s3s2 = -G*self.m2*(self.s3.pos-self.s2.pos)/(pow(self.dist(self.s3.pos,self.s2.pos),3))

    def update(self):
        for i in range(int(1/dt*ddt)): #calculating steps between every visual update
            self.rk4()
        #update positions for visual
        self.s1.pos = self.p1
        self.s2.pos = self.p2
        self.s3.pos = self.p3
        self.forceArrow() #update force arrows

        self.t += ddt #update time

def SetPos(): #initialization of positions and velocity
    global p1,p2,p3,v1,v2,v3
    p1_list = list(map(float,p1_input.text.split(",")))
    p1 = vector(p1_list[0],p1_list[1],p1_list[2])
    p2_list = list(map(float,p2_input.text.split(",")))
    p2 = vector(p2_list[0],p2_list[1],p2_list[2])
    p3_list = list(map(float,p3_input.text.split(",")))
    p3 = vector(p3_list[0],p3_list[1],p3_list[2])
    v1_list = list(map(float,v1_input.text.split(",")))
    v1 = vector(v1_list[0],v1_list[1],v1_list[2])
    v2_list = list(map(float,v2_input.text.split(",")))
    v2 = vector(v2_list[0],v2_list[1],v2_list[2])
    v3_list = list(map(float,v3_input.text.split(",")))
    v3 = vector(v3_list[0],v3_list[1],v3_list[2])

def Set(evt): #drop menu list for stable configuration, setUp
    #fix = IA1
    fix = [[],[]]
    match evt.index:
        case 0:
            fix = [[vector(0,0,0),vector(0,0,0),vector(0,0,0)],[vector(0,0,0),vector(0,0,0),vector(0,0,0)]]
        case 1:
            fix = IA1
        case 2:
            fix = IA2
        case 3:
            fix = IA17
        case 4:
            fix = IA77
        case 5:
            fix = IA115
        case 6:
            fix = IB1
        case 7:
            fix = IB2
        case 8:
            fix = IB6
        case 9:
            fix = IIC1
    m1_input.text = "1.0"
    m2_input.text = "1.0"
    m3_input.text = "1.0"
    p1_input.text = str(fix[0][0].x)+","+str(fix[0][0].y)+","+str(fix[0][0].z)
    p2_input.text = str(fix[0][1].x)+","+str(fix[0][1].y)+","+str(fix[0][1].z)
    p3_input.text = str(fix[0][2].x)+","+str(fix[0][2].y)+","+str(fix[0][2].z)
    v1_input.text = str(fix[1][0].x)+","+str(fix[1][0].y)+","+str(fix[1][0].z)
    v2_input.text = str(fix[1][1].x)+","+str(fix[1][1].y)+","+str(fix[1][1].z)
    v3_input.text = str(fix[1][2].x)+","+str(fix[1][2].y)+","+str(fix[1][2].z)
wtext(text="\nStable configurations: ")
choicelist=['None','IA1','IA2', 'IA17', 'IA77', 'IA115','IB1','IB2','IB6','IIC1'] #choices
menu(bind=Set, choices=choicelist)
wtext(text="\n")

def RanDom(): #create random configuration
    #random mass between 0 and 3
    m1_input.text = random.uniform(0,3)
    m2_input.text = random.uniform(0,3)
    m3_input.text = random.uniform(0,3)
    #random position and velocity between -2 and 2
    v1_input.text = str(random.uniform(-2,2))+","+str(random.uniform(-2,2))+","+str(random.uniform(-2,2))
    v2_input.text = str(random.uniform(-2,2))+","+str(random.uniform(-2,2))+","+str(random.uniform(-2,2))
    v3_input.text = str(random.uniform(-2,2))+","+str(random.uniform(-2,2))+","+str(random.uniform(-2,2))
    p1_input.text = str(random.uniform(-2,2))+","+str(random.uniform(-2,2))+","+str(random.uniform(-2,2))
    p2_input.text = str(random.uniform(-2,2))+","+str(random.uniform(-2,2))+","+str(random.uniform(-2,2))
    p3_input.text = str(random.uniform(-2,2))+","+str(random.uniform(-2,2))+","+str(random.uniform(-2,2))

Random = button(bind=RanDom, text="Random Configuration")
wtext(text="\n\n")

bodies = []
colors = [color.blue, color.red, color.green, color.yellow] #list of colors for every copy
def StartSim(): #Start the simulation
    global p1,p2,p3,v1,v2,v3, bodies
    global Start, Exit, Started
    if not Started:
        SetPos() #define starting position from input
        #define differences between copies from input
        dm_list = list(map(float,dm_input.text.split(",")))
        dp1_list = list(map(float,dp1_input.text.split(",")))
        dp1 = vector(dp1_list[0],dp1_list[1],dp1_list[2])
        dp2_list = list(map(float,dp2_input.text.split(",")))
        dp2 = vector(dp2_list[0],dp2_list[1],dp2_list[2])
        dp3_list = list(map(float,dp3_input.text.split(",")))
        dp3 = vector(dp3_list[0],dp3_list[1],dp3_list[2])
        dv1_list = list(map(float,dv1_input.text.split(",")))
        dv1 = vector(dv1_list[0],dv1_list[1],dv1_list[2])
        dv2_list = list(map(float,dv2_input.text.split(",")))
        dv2 = vector(dv2_list[0],dv2_list[1],dv2_list[2])
        dv3_list = list(map(float,dv3_input.text.split(",")))
        dv3 = vector(dv3_list[0],dv3_list[1],dv3_list[2])
        for i in range(int(nc_input.text)): #create all objects using the inputs given
            bodies.append(ThreeBody(p1+i*dp1,p2+i*dp2,p3+i*dp3,float(m1_input.text)+i*dm_list[0],float(m2_input.text)+i*dm_list[1],float(m3_input.text)+i*dm_list[2], v1+i*dv1,v2+i*dv2,v3+i*dv3, colors[i]))
        Started = True #System has been initialized
    print("Start Simulation")
    Start = True #Begin simulation

def ExitProgram(): #Exit program
    global Start, Exit
    Start = False
    Exit = True
    print("Ending Simulation")

def pause(): #Pause the program
    global Start
    Start = False
    print("Pause...")

StartButton = button(bind=StartSim, text='Start Simulation!') #Button to start the simulation
EndButton = button(bind=ExitProgram, text="Exit!") #Button to exit the simulation
Pause = button(bind=pause, text="Pause") #Button to pause the simulation

def dist(a,b):
    return sqrt(pow(abs(a.s1.pos.x-b.s1.pos.x),2)+pow(abs(a.s1.pos.y-b.s1.pos.y),2)+pow(abs(a.s1.pos.z-b.s1.pos.z),2))

while True:
    rate(1/ddt) #How often it loops per second, limited because of the limited speed at which the computer can render the visual
    if Start == True: #when Start button is pressed
        for b in bodies:
            b.update() #update every copy
        if int(nc_input.text)>=2: #if at least 2 copies, plot graph showing the position difference between the s1 of the first two copies
            axs.append(bodies[0].t)
            ays.append(dist(bodies[0],bodies[1]))
    if Exit:
        break

ax.plot(axs,ays)
ax.set(xlabel="t", ylabel="Distance")
fig.savefig("Distance_TBP.pdf") #save plot as pdf
plt.show() #show the matplotlib plots