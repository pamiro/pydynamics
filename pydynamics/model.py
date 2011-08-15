#    Copyright (c) 2011, PyDynamics Authors
#    All rights reserved.
#    
#    Redistribution and use in source and binary forms, with or without
#    modification, are permitted provided that the following conditions are met:
#        * Redistributions of source code must retain the above copyright
#          notice, this list of conditions and the following disclaimer.
#        * Redistributions in binary form must reproduce the above copyright
#          notice, this list of conditions and the following disclaimer in the
#          documentation and/or other materials provided with the distribution.
#        * Neither the name of the PyDynamics nor the
#          names of its contributors may be used to endorse or promote products
#          derived from this software without specific prior written permission.
#    
#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#    DISCLAIMED. IN NO EVENT SHALL PYDYNAMICS AUTHORS BE LIABLE FOR ANY
#    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
from pydynamics.fwdyn import fwdyn_ab

'''
Created on Jul 29, 2011

@author: Pavel Mironchyk <p.mironchyk at gmail.com>
'''

import numpy    as np
import spatvect as sv
from copy import copy

class RigidBody:
    """Rigid body with CoM and spacial Inertia """
    def __init__(self, mass, CoM, mI):         
        #joints that connect body to other bodies"""
        self.joints  = {}
        self.mass = mass
        self.CoM  = CoM        
        # spacial rigid body inertia
        self.rbI  = sv.mcI(mass, CoM, mI)        
        # moment of inertia tensor
        self.mI   = mI
        # External force
        self.f_ext = None

        self.figure = None
        
        self.Xup = None
        
        self.descr = "none"
        
    def position(self):
        return sv.vector3d(self.Xup) + np.matrix(self.CoM)
        
    def description(self):
        return self.descr
    
    def describe(self, text):
        self.descr = text  
      
class Joint:
    def __init__(self, transform, a, b):
        self.bodies = [a, b]
        a.joints[b] = self  
                      
        # Join Transformation
        self.Xj = None
        # Joint Space
        self.S = None    
        # Frame A to frame B transformation 
        self.Xtree = transform
        
        self.figure = None
        
        self.tau = 0.0

    def integrate_q(self, qd, step):
        self.q += qd*step

    def description(self):
        return self.descr
    
    def describe(self, text):
        self.descr = text


class GenericModel:
    def __init__(self):
        self.bodies = []
        self.joints = []
        
        # model gravity
        self.a_grav = np.matrix([0,0,0,0,0, 9.81]).transpose(); 
 
    def addbody(self, b):
        self.bodies.append(b)
    
    def addjoint(self, j):
        for body in j.bodies:
            if not (body in self.bodies):
                self.bodies.append(body)
        self.joints.append(j)
    
    def add(self, o):
        if isinstance(o, RigidBody):
            self.addbody(o)
        if isinstance(o, Joint):
            self.addjoint(o)

class KTree(GenericModel):
    """Kinematic tree""" 
     
    def __init__(self):
        GenericModel.__init__(self)
        
    def update_parent(self):
        self.parent = []
        if len(self.bodies) > 0:
            for bindex in range(1, len(self.bodies)):
                body = self.bodies[bindex]
                for joint in self.joints:
                    if body in joint.bodies:        
                        ind = (joint.bodies.index(body) + 1) % 2
                        other_body = joint.bodies[ind]
                        other_body_ind = self.bodies.index(other_body)        
                        if other_body_ind < bindex:
                            self.parent.append(other_body_ind)
        return self.parent
        
 
    def update_representation(self):
        self.predecessor = []
        self.successor = []
        for j in self.joints:
            self.predecessor.append(self.bodies.index(j.bodies[0])) 
            self.successor.append(self.bodies.index(j.bodies[1]))    
            
        return (self.predecessor, self.successor)

    def update(self):
        self.update_parent()
        self.update_representation


class KClosedLoop(GenericModel):
    """Kinematic closed loop"""
    pass


class World:
    def __init__(self):
        self.models = []
        self.base = RigidBody(0.0, [0.0, 0.0, 0.0], np.matrix([[0, 0, 0],
                                                                [0, 0, 0],
                                                                [0, 0, 0]]))
        self.base.Xup = sv.Xtrans([0.0, 0.0, 0.0])

        self.base.describe("world base")

    def add(self, model):        
        self.models.append(model)
        
    def update(self):
        for model in self.models: model.update()
        
    def simulate(self, step):
        for model in self.models: 
            fwdyn_ab(model)
            for joint in model.joints:
                joint.qd += joint.qdd * step
                joint.integrate_q(joint.qd, step)