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

import math
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


def rbda_eq_4_12(qe):
    p0, p1, p2, p3 = qe[0,0], qe[1,0], qe[2,0], qe[3,0]
    return np.matrix([[p0**2+p1**2-0.5,  p1*p2+p0*p3,     p1*p3-p0*p2,],
                      [p1*p2-p0*p3,      p0**2+p2**2-0.5, p2*p3+p0*p1 ],
                      [p1*p3+p0*p2,      p2*p3-p0*p1,     p0**2+p3**2-0.5]]) * 2.0 

def rbda_eq_4_13(qe):
    p0, p1, p2, p3 = qe[0,0], qe[1,0], qe[2,0], qe[3,0]
    return np.matrix([[-p1,  -p2, -p3],
                      [ p0,  -p3,  p2],
                      [ p3,   p0, -p1],
                      [ -p2,  p1,  p0]]) * 0.5

def normalize(qe):
    sum = 0.0; 
    for i in xrange(4): sum += qe[i,0]*qe[i,0]
    qe_norm = qe / math.sqrt(sum)
    return qe_norm

class DoF0Joint(Joint):
    def __init__(self, bodyA, bodyB):
        Joint.__init__(self, bodyA, bodyB)
        
        self.Xj = sv.Xrot(np.matrix([[1.0,0.0,0.0],
                                     [0.0,1.0,0.0],
                                     [0.0,0.0,1.0]]))*sv.Xtrans([0.0,0.0,0.0])
        self.S  = np.matrix([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]).transpose()

    def update(self, q):
        pass
    
        
class DoF6Joint(Joint):
    def __init__(self, transform, a, b, qe = None, qr = None):
        Joint.__init__(self, transform, a, b)
        self.qe = qe
        self.qr = qr
        if qe is None:
            theta = 0.0
            ux    = 0.0
            uy    = 0.0
            uz    = 0.0
            p0 = np.cos(theta/2)
            p1 = np.sin(theta/2)*ux
            p2 = np.sin(theta/2)*uy
            p3 = np.sin(theta/2)*uz 
            self.qe = np.matrix([p0, p1, p2, p3]).transpose()
        if qr is None:
            self.qr = np.matrix([0.0,0.0,0.0]).transpose()
        
        self.qe_norm = normalize(self.qe)
        self.qd =  np.resize(0.0, (6,1))
        self.qdd = np.resize(0.0, (6,1))
        self.describe(" 6DoF Joint")
        self.update(self.qe, self.qr)
        
    def update(self, qe, qr):
        self.qe = qe
        self.qr = qr
        self.e = rbda_eq_4_12(qe)
        self.r = -np.linalg.inv(self.e)*qr
        self.Xj = sv.Xrot(self.e)*sv.Xtrans(self.r)
        self.S  = sv.mkdiagm(np.ones(6))

    # thanks scbtx authors for the hint
    def integrate_q(self, qd, dt):
        w_body_frame, v_body_frame = (qd[0:3], qd[3:])
        qed = rbda_eq_4_13(self.qe_norm) * w_body_frame
        new_qe = normalize(self.qe + qed * dt)
        qrd = self.e.transpose() * v_body_frame
        new_qr = self.qr + qrd * dt
        self.update(new_qe, new_qr)
    
    