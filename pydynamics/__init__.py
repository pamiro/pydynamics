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
This object-oriented library for rigid body dynamics simulations. 

Modules include:

    :mod:`pydynamics.sv`
        defines various routines of spatial vector algebra

    :mod:`pydynamics.model`
        defines the :class:`~pydynamics.model.KTree` class.

    :mod:`pydynamics.joints`
        defines various types of joints that can be used in models

    :mod:`pydynamics.minert`
        defines various functions for constructing moment of intertia tensors

    :func:`pydynamics.invdyn`
        define function for calculating inverse dynamics of a model
    
    :func:`pydynamics.fwdyn_ab`
        define function for calculating forward dynamics of a model (kinematic tree, articulated body)
    

Created on Jul 29, 2011

@author: Pavel Mironchyk <p.mironchyk at gmail.com>
'''

import spatvect      as sv
import model         
import joints        
import visualization as visual
from pydynamics.inertia import InertiaTensor as minert
from fwdyn  import fwdyn_ab
from invdyn import invdyn
import numpy as np
