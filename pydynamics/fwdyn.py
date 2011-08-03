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
Created on Aug 3, 2011

@author: Pavel Mironchyk <p.mironchyk at gmail.com>
'''

import numpy    as np
import spatvect as sv

def fwdyn_ab(model):
    """Forward dynamics of articulated body"""

    qdd = []    
    for i in range(0, len(model.parent)):
        pbody_ind = model.parent[i]
        pbody = model.bodies[pbody_ind]
        body  = model.bodies[i + 1]
        joint = pbody.joints[body]
                   
        vJ = joint.S*joint.qd;
        joint.Xup = joint.Xj * joint.Xtree;
        
        if model.parent[i] == 0:
            body.v = vJ;
            body.c = np.zeros((6,1))
        else:
            body.v = joint.Xup*pbody.v + vJ;
            body.c = sv.cross_m_m(body.v, vJ)

        body.IA = body.rbI
        body.pA = sv.cross_m_f(body.v, body.rbI*body.v) 
        
        if body.f_ext != None:
            body.pA = body.pA - body.f_ext

            
    for i in range(len(model.parent)-1,-1, -1):
        pbody_ind = model.parent[i]
        pbody = model.bodies[pbody_ind]
        body  = model.bodies[i + 1]
        joint = pbody.joints[body]

        joint.U = body.IA * joint.S
        joint.d = joint.S.transpose() * joint.U
        joint.u = joint.tau - (joint.S.transpose() * body.pA)

        if model.parent[i] != 0:
            Ia = body.IA - (joint.U/joint.d*joint.U.transpose())
            pa = body.pA + (Ia*body.c) + (joint.U * joint.u / joint.d)
            pbody.IA = pbody.IA + (joint.Xup.transpose() * Ia * joint.Xup)
            pbody.pA = pbody.pA + (joint.Xup.transpose() * pa)
            

    for i in range(0, len(model.parent)):

        pbody_ind = model.parent[i]
        pbody = model.bodies[pbody_ind]
        body  = model.bodies[i + 1]
        joint = pbody.joints[body]

        if model.parent[i] == 0:
            body.a = joint.Xup * -model.a_grav + body.c
        else:
            body.a = joint.Xup * pbody.a + body.c
        joint.qdd = (joint.u - joint.U.transpose()*body.a)/joint.d;
        body.a = body.a + joint.S * joint.qdd
        qdd.append(joint.qdd)
        
    return qdd
