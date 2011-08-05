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
Routines to support spacial vectors algebra

Created on May 22, 2011

@author: Pavel Mironchyk <p.mironchyk at gmail.com>
'''

import numpy as np

def Xrotx(theta):
    """spacial rotational transformation around X axes for angle theta"""
    c = np.cos(theta);
    s = np.sin(theta);
    return np.matrix([[1, 0, 0, 0, 0, 0 ],
                      [0, c, s, 0, 0, 0 ],
                      [0,-s, c, 0, 0, 0 ],
                      [0, 0, 0, 1, 0, 0 ],
                      [0, 0, 0, 0, c, s ],
                      [0, 0, 0, 0,-s, c]]);


def Xroty(theta):
    """spacial rotational transformation around Y axes for angle theta"""
    c = np.cos(theta);
    s = np.sin(theta);
    return np.matrix([[c, 0,-s, 0, 0, 0],
                      [0, 1, 0, 0, 0, 0],
                      [s, 0, c, 0, 0, 0],
                      [0, 0, 0, c, 0,-s],
                      [0, 0, 0, 0, 1, 0],
                      [0, 0, 0, s, 0, c]]);


def Xrotz(theta):
    """spacial rotational transformation around Z axes for angle theta"""
    c = np.cos(theta);
    s = np.sin(theta);
    return np.matrix([[c, s, 0, 0, 0, 0],
                     [-s, c, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, c, s, 0],
                      [0, 0, 0,-s, c, 0],
                      [0, 0, 0, 0, 0, 1]])

def Xrot(E):
    """spacial rotational transformation defined by rotational matrix E"""
    m = np.resize(0.0, (6,6))
    m[0:3,0:3] = E
    m[3:6,3:6] = E
    return m

def Xtrans(r):
    """spacial linear transformation defined by 3D displacement vector r"""
    return np.matrix([
                      [ 1,     0,     0,     0,   0,   0 ],
                      [ 0,     1,     0,     0,   0,   0 ],
                      [ 0,     0,     1,     0,   0,   0 ],
                      [ 0,     r[2], -r[1],  1,   0,   0 ],
                      [-r[2],  0,     r[0],  0,   1,   0 ],
                      [ r[1], -r[0],  0,     0,   0,   1 ]])


def cross_m_m(v, m):
    """cross product between two motion vectors""" 
    result = np.resize(0.0, (6,1))
    result[0:3,0] = np.cross(v[0:3].transpose(), m[0:3].transpose())
    result[3:6,0] = np.cross(v[0:3].transpose(), m[3:6].transpose()) + np.cross(v[3:6].transpose(), m[0:3].transpose())
    return result


def cross_m_f(v, f):
    """cross product between motion and force vectors"""
    result = np.resize(0.0, (6,1))
    result[0:3,0] = np.cross(v[0:3].transpose(), f[0:3].transpose()) +  np.cross(v[3:6,0].transpose(), f[3:6,0].transpose())
    result[3:6,0] = np.cross(v[0:3].transpose(), f[3:6].transpose())
    return result


def mcI(mass, CoM, I):
    """spacial rigid body inertia from mass, CoM and moment of inertia tensor"""
    c = CoM
    C = np.matrix([[ 0,    -c[2],  c[1]],
                   [ c[2],  0,    -c[0]],
                   [-c[1],  c[0],  0 ]])

    rbi = np.resize( 0.0 ,(6,6))
    rbi[0:3,0:3] = I + mass*C*C.transpose()
    rbi[0:3,3:6] = mass*C
    rbi[3:6,0:3] = mass*C.transpose()
    rbi[3:6,3:6] = mass*np.eye(3)
    return rbi

def mkdiagm(arr):
    nmarr = np.array([])
    nmarr.resize((len(arr), len(arr)))
    for i in range(0, len(arr)):
        nmarr[i,i] = arr[i]  
    return np.matrix(nmarr)

def vector3d(mvector):
    """vector representation of spacial transform"""
    rcrosst = mvector[0:3,0:3].transpose() * mvector[3:6,0:3];
    rx = rcrosst[1,2];
    ry = rcrosst[2,0];
    rz = rcrosst[0,1];    
    return [rx, ry, rz]