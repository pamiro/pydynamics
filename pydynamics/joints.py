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

'''
Created on Jul 29, 2011

@author: Pavel Mironchyk <p.mironchyk at gmail.com>
'''

import numpy    as np
import spatvect as sv
from model import Joint

class RevoluteJoint(Joint):
    def __init__(self, bodyA, bodyB, q):
        Joint.__init__(self, bodyA, bodyB)
        self.update(q)
    
    def update(self, q):
        self.Xj = sv.Xrotz(q)
        self.S = np.matrix([0,0,1,0,0,0]).transpose()

class PrismaticsJoint(Joint):
    def __init__(self, bodyA, bodyB, d):
        Joint.__init__(self, bodyA, bodyB)
        self.update(d)
    
    def update(self, d):
        self.Xj = sv.Xtrans(np.array([0, 0, d]))
        self.S  = np.matrix([0,0,0,0,0,1]).transpose();

class HelicalJoint(Joint):
    def __init__(self, bodyA, bodyB, q, d):
        Joint.__init__(self, bodyA, bodyB)
        self.update(q, d)
    
    def update(self, q, d):
        self.Xj = sv.Xrotz(q) * sv.Xtrans(np.array([0, 0, q*d]));
        self.S  = np.matrix([0, 0, 1, 0, 0, d]).transpose();

class DoF6Joint(Joint):
    def __init__(self, bodyA, bodyB, q):
        Joint.__init__(self, bodyA, bodyB)
        if q == None:
            self.q = np.zeros((6,1))
        else:
            self.q = q
        self.qd = np.zeros((6,1))
        self.qdd = np.zeros((6,1))
        self.update(self.q)
        
    def update(self, q):
        e = sv.Xrotx(q[2,0])*sv.Xroty(q[1,0])*sv.Xrotz(q[0,0])
        E = (e)[0:3,0:3]
        r = -np.linalg.inv(E)*(np.matrix([q[3,0],q[4,0],q[5,0]]).transpose())
        self.Xj = e*sv.Xtrans(r)
        self.S  = sv.mkdiagm(np.ones(6))
