'''
Created on Aug 2, 2011

@author: Pavel Mironchyk <p.mironchyk at gmail.com>
'''

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot        as plt
import numpy                    as np
import pydynamics.model         as m
import pydynamics.joints        as j
import pydynamics.visualization as v
import pydynamics.spatvect      as sv
from   pydynamics.invdyn import invdyn
from   pydynamics.fwdyn  import fwdyn_ab
import math
import os

model = m.KTree()

def mkdiagm(arr):
    nmarr = np.array([])
    nmarr.resize((len(arr), len(arr)))
    for i in range(0, len(arr)):
        nmarr[i,i] = arr[i]  
    return np.matrix(nmarr)

inf = 1.0e+99

q = [ 1.97746724216702, 
      2.54966593680547,
      -2.34371095524892,
      2.59731710563547,
      0.831637671152858,
      -2.52872821404116,
      -1.39173673672867,
      0.294565272647012,
      2.87460022633501,
      2.92098081384033,
      -2.15128045457413]

qd = [ 0.914333896485891,
      -0.0292487025543176,
      0.600560937777600,
      -0.716227322745569,
      -0.156477434747450,
      0.831471050378134,
      0.584414659119109,
      0.918984852785806,
      0.311481398313174,
      -0.928576642851621,
      0.698258611737554]

qdd = [0.357470309715547,
       0.515480261156667,
       0.486264936249832,
       -0.215545960931664,
       0.310955780355113,
       -0.657626624376877,
       0.412092176039218,
       -0.936334307245159,
       -0.446154030078220,
       -0.907657218737692,
       -0.805736437528305]

nb    = 12
bf    = 1.0
skew  = 45.0
taper = 0.95

length = []
parent = []

pitch = np.zeros((nb,1))
pitch[3 - 1] = 0.1;
pitch[5 - 1] = inf;
pitch[7 - 1] = -0.1;
pitch[9 - 1] = inf;

joints = []

for i in range(0,nb):                            
    parent.append(int(math.floor((i-1+math.ceil(bf))/bf)))
    length.append(taper**(i));
    
    body = m.RigidBody()
    body.CoM = length[i] * np.array([0.5,0,0]);
    body.mass = taper**(3*(i));
    Icm = body.mass * (length[i]**2) * mkdiagm([0.0025,1.015/12,1.015/12]);
    body.rbI = sv.mcI(body.mass, body.CoM, Icm)
    body.ii = i
    model.addbody(body)        

for i in range(0,len(parent)-1):
    
    if parent[i] == 0:
        Xtree = sv.Xtrans([0, 0, 0]);
    else:
        l = length[parent[i]-1]
        Xtree = sv.Xrotx(skew)*sv.Xtrans([l,0,0])
        
    joint = None
    if pitch[i] == 0:
        joint = j.RevoluteJoint(model.bodies[parent[i]],model.bodies[i+1], q[i])
    else:
        if pitch[i] == inf: #% prismatic joint
            joint = j.PrismaticsJoint(model.bodies[parent[i]],model.bodies[i+1], q[i])
        else: 
            joint = j.HelicalJoint(model.bodies[parent[i]],model.bodies[i+1], q[i], pitch[i])
    joint.Xtree = Xtree
    joint.qd = qd[i]
    joint.qdd = qdd[i]
    if i < 11:
        model.parent = model.parent[0:11] 
        model.add(joint)
       
model.update()

tau = invdyn(model)
new_qdd = fwdyn_ab(model)

for i in range(0, len(new_qdd)):
    if math.fabs(qdd[i] - new_qdd[i]) > 1e-13:
        print "Test failed"
        os._exit(-1) 
    else:
        print qdd[i]- new_qdd[i]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_autoscale_on(False)
ax.set_xlim3d(-3,3);
ax.set_ylim3d(-3,3);
ax.set_zlim3d(-3,3);
#ax.cla()
v.plotktree(ax,model)
plt.show()
