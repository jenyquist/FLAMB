# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
from __future__ import division
from collections import namedtuple
from math import *
from scipy.special import expn
import numpy as np

# ------------------------------------------------------------------------------
def parameters():
    '''This function is the place to set all of the parameters.'''
    
    #PARAMETER DEFINITIONS
    
    # 1.1 Borehole parameters
    
    # Distance between boreholes [meters]
    # (approximated as similar in each aquifer)
    L = 21
    
    # Radius of the pumped (i=0) and observation (i=1) boreholes [meters]
    R = [3.75e-2, 3.75e-2]
    
    # 1.2 Pumping Parameters
    
    # Pumping time [minutes * 60 to give seconds]
    pump_time = 20 * 60
    
    # Pumping flow rate [L/min]
    Q = 17.8
    
    # 1.3 Definition of the aquifers
    # Number of aquifers
    Naq = 5
    
    # Depth in [m] from surface of each aquifer in pumping borehole
    depth0 = [26, 56, 96, 116, 142]
    
    # Depth in [m] from the surface of each aquifer in the observation borehole
    depth1 = [26, 52, 92, 116, 139]
    
    # Transmissivities of the aquifers [m2/sec]
    Trans_aq = 1/(24*3600)*np.array([1.8, 2, 2, 0.2, 0.2])
    
    # Storativities of the aquifers [dimensionless]
    Stora_aq = np.array([1e-5, 1e-5, 1e-5, 1e-5, 1e-5])
    
    # Aquifer connections
    pump_well_connections = [1]
    obs_well_connections = [0, 1, 2, 3, 4]
    Aq_well_connect = [pump_well_connections, obs_well_connections]
                       
    # 1.4 Time discretization
    
    # Maximum time of the simulation [sec]
    t_max = 40*60
    
    # Number of time steps
    Nt = 100
    
    # 2. MODEL PARAMETER COMPUTATION
    
    # 2.1 Time Discretization
    delta_t = t_max/Nt
    vect_t = np.arange(1,Nt+1)*delta_t
    
    # 2.2 Distances between aquifers
    
    # Pumped Borehole
    dist_aq0 = np.array(depth0)
    dist_aq0[1:] = dist_aq0[1:]-dist_aq0[0:-1]

    
    # Observation borehole
    dist_aq1 = np.array(depth1)
    dist_aq1[1:] = dist_aq1[1:]-dist_aq1[0:-1]
    
    dist_aq = np.array([dist_aq0, dist_aq1])

    # Flow velocity at the top of the pumped well [m/s]
    q_pumped = Q*1e-3/(60*np.pi*R[0]**2)
    
    params = namedtuple('params',['L', 'R', 'pump_time', 'Trans_aq', 'Stora_aq',
                    'vect_t', 'dist_aq', 'q_pumped', 'Naq', 'Aq_well_connect'])
    p = params(L,R,pump_time,Trans_aq,Stora_aq,vect_t,dist_aq,q_pumped,Naq,Aq_well_connect)
    return p

#-------------------------------------------------------------------------------
def return_conductivities(R,dist_aq):
    '''Returns conductivities normalized by the distance between aquifers'''
    # Density of water [kg/m3]
    RHO = 1.0e6 
    # Acceleration of gravity on the Earth's surface [m/s2]
    G = 9.81
    # Viscosity of water under standard temperature and pressure [MPa-sec]
    MU = 1.0
    # Conductivity normalized by the distance between aquifers
    betaI = np.zeros(dist_aq.shape)
    betaI[0,:] = RHO*G*R[0]**2/(8*MU*dist_aq[0,:]) # pumped well
    betaI[1,:] = RHO*G*R[1]**2/(8*MU*dist_aq[1,:]) # observation well
    return betaI

# ------------------------------------------------------------------------------
    
def return_flow_indices(Aq_well_connect, Nt, Naq):
    ''' Return the number for the flow in the linear system for a given 
    well i and aquifer num_aq'''
        
    Indices = np.zeros((len(Aq_well_connect), Naq, Nt))
    cpt_flow = 1  # Number of connections
    
    # Loop over the wells
    for i in np.arange(len(Aq_well_connect)):
        # Loop over the aquifers
        for I in np.arange(len(Aq_well_connect[i])):
            num_aqi = Aq_well_connect[i][I]
            # Loop over time
            for k in np.arange(Nt):
                Indices[i][num_aqi][k] = cpt_flow
                cpt_flow += 1
    cpt_flow -= 1
        
    return (Indices, cpt_flow)
    
#------------------------------------------------------------------------------

def return_boreholes_list(Aq_well_connect, Naq):
    ''' Returns the list of aquifers with boreholes intersecting for aquifer'''  
    WellsInterAq = [[] for i in range(Naq)]
    for i in range(Naq):
        for j in range(len(Aq_well_connect)):
            if (i in Aq_well_connect[j]):
               WellsInterAq[i].append(j)
    return WellsInterAq 

# ------------------------------------------------------------------------------

def flow_velocity_computation(p):
    ''' Define the linear system of equations AX=b where X = [X0 X1] with X0
    and X1 the flow velocity in the pumped and observation borehole,
    respectively.'''
    
    # 0.0 PARAMETERS AND INPUTS
    
    # 0.1 Aquifer parameters
    Nt = len(p.vect_t)
    delta_t = p.vect_t[1] - p.vect_t[0]
    
    # 0.2 Physical parameters
    
    # Density of water [kg/m3]
    rho = 1.0e6 
    
    # Acceleration of gravity on the Earth's surface [m/s2]
    g = 9.81
    
    # Viscosity of water under standard temperature and pressure [MPa-sec]
    mu = 1.0
    
    # 0.3 Conductivity normalized by the distance between aquifers
    betaI = return_conductivities(p)
    
    # 1.0 LINEAR SYSTEM DEFINITION
    
    # Flow indices
    Indices, cpt_flow = return_flow_indices(p.Aq_well_connect, Nt, p.Naq)
    
    # List of boreholes that intersect each well
    WellsInterAq = return_boreholes_list(p.Aq_well_connect, p.Naq)
    
    # Computation of the integrals required in the linear system
    
    # WORKING ON THIS
    IntijI = IntergralComputation(p)
    
    return (WellsInterAq)
    
# ------------------------------------------------------------------------------
    
def fun_G1(r,theta,t1,t2,xij,alpha):
    beta = ((xij-r*cos(theta))**2 + (r*sin(theta))**2)/(4*alpha)
    if t1==0:
        res = r*expn(1,beta/t2)
    else:
        res = r*(expn(1,beta/t2) - expn(1,beta/t1))      
    return res

# ------------------------------------------------------------------------------

def fun_G2(r,t,alpha):
    c = 4 * alpha * t
    if t==0:
        res = 0.0
    else:
        res = c*pi*(1 - exp(-r**2/c)) + pi*r**2*expn(1,r**2/c)
    return res

# ------------------------------------------------------------------------------
# IntijI=IntegralComputation(Aq_well_connect,vect_t,L,Trans_aq,Stora_aq,R);
def intergral_computation(p):
    ''' Function to compute the integrals required in the linear system'''
    
    # Variables
    nwell = len(p.Aq_well_connect)
    alpha = p.Trans_aq/p.Stora_aq # Ratio be transmisivity and storativity
    delta_t = p.vect_t[1] - p.vect_t[0] # Timestep
    
    # Arrays to store the computed integrals for the wells i,j and 
    # aquifer I (must pass shape to zeros as a tuple)
    IntijI = np.zeros((nwell, nwell, p.Naq, len(p.vect_t)))
    
    # J.E. NYQUIST: Note the original matlab code was written here to handle
    # more than two wells, but never tested for that.  I've simplified by
    # hardwiring the two well case.  This will need to be changed to add
    # additional wells.
    # For each well
    
    for i in range(len(nwell)):
        for I in range(len(p.Aq_well_connect[0])):
            num_aq = p.Aq_well_connect[0][I]
            # Second loop on the wells to check if the aquifer intersects both
            for j in range(len(nwell)):
                if (num_aq in p.Aq_well_connect[1]):
                    pass
    return IntijI

if __name__ == "__main__":
    # 
    #  MAIN ROUTINE
    #                   
    # ------------------------------------------------------------------------------
    # -*- coding: utf-8 -*-
    '''This program is a Python translation of Matlab routines written by Delphine
    Roubinet.  The methodology is documented in:
        
    Roubinet, D., Irving, J., & Day-lewis, F. D. (2015). Advances in Water Resources
    Development of a new semi-analytical model for cross-borehole flow experiments 
    in fractured media. Advances in Water Resources, 76, 97–108. 
    doi:10.1016/j.advwatres.2014.12.002
    
    Python Code by Jonathan E. Nyquist
    February, 2015
    '''
    
    # Load the module with all the functions
    #import xflow
    
    # 1. PARAMETER DEFINITION
    # All of the model parameters are set in this function
    # p is a named tuple with the following parameters:
    # L,R,pump_time,Trans_aq,Stora_aq,vect_t,dist_aq,q_pumped,Naq,Aq_well_connect
    # Each can be accessed with the dot notation (eg. p.pump_time)
    p = parameters()
    
    # 2. FLOW VELOCITY DETERMINATION
    # 2.1 Flow velocity computation from the developed forward model
    # q_pump, q_obs = flow_velocity_computation(p)
    WellsInterAq = flow_velocity_computation(p)
    print(WellInterAq)
    #
    
    # 2.2 Conversion from [m/s] to [L/min]
    
    # 3. PLOT RESULTS FOR FLOW VELOCITY COMPARISON
