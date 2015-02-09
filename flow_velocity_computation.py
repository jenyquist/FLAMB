from __future__ import division
import numpy as np

def flow_velocity_computation(p):
    ''' Define the linear system of equations AX=b where X = [X0 X1] with X0
    and X1 the flow velocity in the pumped and observation borehole,
    respectively.'''
    
    # 0. PARAMETERS AND INPUTS
    
    # 0.1 Aquifer parameters
    Nt = len(p.vect_t)
    delta_t = p.vect_t[1] - p.vect_t[0]
    
    return delta_t