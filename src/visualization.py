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
Few routines to visualize simple geometric figures using matplotlib

Created on Jul 21, 2011

@author: Pavel Mironchyk Pavel Mironchyk <pmironchyk at gmail.com>
'''

from copy import copy
import numpy as np

class CubeFigure:
    def __init__(self, w, h, d):
        u = np.linspace(-w/2, w/2, 4)
        v = np.linspace(-h/2, h/2, 4)
        p = np.linspace(-d/2, d/2, 4)
        
        self.x,self.y = np.meshgrid(u,v)
        self.z = d/2 * np.outer(np.ones(np.size(u)), np.ones(np.size(v)))
        self.z0 = self.z - d;
        
        self.x1,self.z1 = np.meshgrid(u,p)
        self.y1 = h/2 * np.outer(np.ones(np.size(u)), np.ones(np.size(p)))
        self.y0 = self.y1 - h;

        self.y2,self.z2 = np.meshgrid(v,p)
        self.x2 = w/2 * np.outer(np.ones(np.size(u)), np.ones(np.size(p)))
        self.x0 = self.x2 - w;
        
    
    def draw(self, ax, transform = None):
        if transform == None:
            def transform(x,y,z): return (x,y,z);         
        
        x=self.x;
        y=self.y;
        z=self.z;
        z0=self.z0
         
        xt = copy(x)
        yt = copy(y)
        zt = copy(z)
    
        
        for i  in range(0,len(x)):
            for j  in range(0,len(x[i])):
                xt[i][j], yt[i][j], zt[i][j] = transform(x[i][j], y[i][j], z[i][j])
        ax.plot_surface(xt, yt, zt, color='g')    

        for i  in range(0,len(x)):
            for j  in range(0,len(x[i])):
                xt[i][j], yt[i][j], zt[i][j] = transform(x[i][j], y[i][j], z0[i][j])          
        ax.plot_surface(xt, yt, zt, color='g')

        x=self.x1;
        y=self.y1;
        z=self.z1;
        y0=self.y0;
         
        xt = copy(x)
        yt = copy(y)
        zt = copy(z)
    
        for i in range(0,len(x)):
            for j  in range(0,len(x[i])):
                xt[i][j], yt[i][j], zt[i][j] = transform(x[i][j], y[i][j], z[i][j])
        ax.plot_surface(xt, yt, zt, color='r')

        for i  in range(0,len(x)):
            for j  in range(0,len(x[i])):
                xt[i][j], yt[i][j], zt[i][j] = transform(x[i][j], y0[i][j], z[i][j])
        ax.plot_surface(xt, yt, zt, color='r')

        x=self.x2;
        y=self.y2;
        z=self.z2;
        x0=self.x0;
         
        xt = copy(x)
        yt = copy(y)
        zt = copy(z)

        for i in range(0,len(x)):
            for j  in range(0,len(x[i])):
                xt[i][j], yt[i][j], zt[i][j] = transform(x[i][j], y[i][j], z[i][j])
        ax.plot_surface(xt, yt, zt, color='b')

        for i  in range(0,len(x)):
            for j  in range(0,len(x[i])):
                xt[i][j], yt[i][j], zt[i][j] = transform(x0[i][j], y[i][j], z[i][j])
        ax.plot_surface(xt, yt, zt, color='b')


class CylinderFigure:
    def __init__(self, r, h):
        u = np.linspace(0, 2*np.pi, 40)
        x = r*np.sin(u)
        y = r*np.cos(u)
        self.X,A = np.meshgrid(x,y)
        self.Y,B = np.meshgrid(y,x)
        self.Z = np.outer(np.linspace(-h/2,h/2,40), np.ones(np.size(u)))

    def draw(self, ax, transform = None):
        xt = copy(self.X)
        yt = copy(self.Y)
        zt = copy(self.Z)
        
        if transform == None:
            def transform(x,y,z): return (x,y,z);         
        
        for i  in range(0,len(self.X)):
            for j  in range(0,len(self.X[i])):
                xt[i][j], yt[i][j], zt[i][j] = transform(self.X[i][j], 
                                                         self.Y[i][j], 
                                                         self.Z[i][j])
        ax.plot_surface(xt, yt, zt,  color='b')


class Ellipsoid:
    def __init__(self, a, b, c):
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        self.x = a * np.outer(np.cos(u), np.sin(v))
        self.y = b * np.outer(np.sin(u), np.sin(v))
        self.z = c * np.outer(np.ones(np.size(u)), np.cos(v))

    def draw(self, ax, transform = None):
        x=self.x;
        y=self.y;
        z=self.z;
        xt = copy(x)
        yt = copy(y)
        zt = copy(z)
    
        if transform == None:
            def transform(x,y,z): return (x,y,z);         
        
        for i  in range(0,len(x)):
            for j  in range(0,len(x[i])):
                xt[i][j], yt[i][j], zt[i][j] = transform(x[i][j], y[i][j], z[i][j])

        ax.plot_surface(xt, yt, zt, color='r')    


class Sphere(Ellipsoid):
    def __init__(self, r):
        Ellipsoid.__init__(self, r, r, r)

