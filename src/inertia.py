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
Created on Jul 20, 2011

@author: Pavel Mironchyk <pmironchyk at gmail.com>
'''

import numpy as np

class InertiaTensor:
    """Definitions of various inertia tensors,
taken from http://en.wikipedia.org/wiki/List_of_moment_of_inertia_tensors"""
    
    @staticmethod
    def solid_sphere(r, m):
        """inertia tensor for solid sphere of radius r and mass m"""
        e = 2.0/5.0*m*(r**2)
        return np.matrix([[e, 0, 0],
                          [0, e, 0],
                          [0, 0, e]])
    
    @staticmethod
    def solid_ellipsoid(a,b,c, m):
        """inertia tensor for solid ellipsoid of semi-axes a,b,c and mass m"""
        return np.matrix([[1.0/5.0*m*((b**2)+(c**2)), 0, 0],
                          [0,  1.0/5.0*m*((a**2)+(c**2)), 0],
                          [0, 0, 1.0/5.0*m*((b**2)+(a**2))]])
    
    @staticmethod
    def solid_cuboloid(w, h, d, m):
        """inertia tensor for solid cube(loid) of weight w, height h, depth d and mass m"""
        return np.matrix([[1.0/12.0*m*((h**2)+(d**2)), 0, 0],
                          [0,  1.0/12.0*m*((w**2)+(d**2)), 0],
                          [0, 0, 1.0/12.0*m*((w**2)+(h**2))]])
    
    @staticmethod
    def solid_cylinder(r, h, m):
        """inertial tensor for solid cylinder of radius r, height h and mass m"""
        return np.matrix([[1.0/12.0*m*((3*r**2)+(h**2)), 0,  0],
                          [0,  1.0/12.0*m*((3*r**2)+(h**2)), 0],
                          [0,  0,             1.0/2.0*m*(r**2)]])

