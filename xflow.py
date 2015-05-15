# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
from __future__ import division
from collections import namedtuple
from math import *
import numpy as np
from scipy.special import exp1
from scipy.integrate import dblquad
import pdb

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
    pump_time = 2 * 60
    
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
    # WILL NEED TO CHANGE THESE BACK
    t_max = 2*60
    
    # Number of time steps
    Nt = 2
    
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
def return_distance_wells(i, j, L):
    if i==j:
        xij = 0
    elif i==0 and j==1:
        xij = -L
    elif i==1 and j==0:
        xij = L
    else:
        print('WARNING IN RETURN_DISTANCE_WELLS: distance not defined')
        xij = -1
    return xij
# ------------------------------------------------------------------------------
def return_flow_indices(Aq_well_connect, Nt, Naq):
    ''' Return the number for the flow in the linear system for a given 
    well i and aquifer num_aq'''
        
    Indices = np.zeros((len(Aq_well_connect), Naq, Nt))
    cpt_flow = 0  # Number of connections Matlab starts at 1
    
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

def indices(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]

# ------------------------------------------------------------------------------
# function [q_pump,q_obs]=FlowVelocityComputation(vect_t,Trans_aq,Stora_aq,dist_aq,R,L,q_pumped,Aq_well_connect,pump_time,Naq)
def flow_velocity_computation(vect_t,Trans_aq,Stora_aq,dist_aq,R,L,q_pumped,Aq_well_connect,pump_time,Naq):
    import pdb
    ''' Define the linear system of equations AX=b where X = [X0 X1] with X0
    and X1 the flow velocity in the pumped and observation borehole,
    respectively.'''
    
    # 0.0 PARAMETERS AND INPUTS
    
    # 0.1 Aquifer parameters
    Nt = len(vect_t)
    delta_t = vect_t[1] - vect_t[0]
    
    # 0.2 Physical parameters are defined in return_conductivities    
    # 0.3 Conductivity normalized by the distance between aquifers
    betaI = return_conductivities(R, dist_aq)
 
    # 1.0 LINEAR SYSTEM DEFINITION
    
    # Flow indices
    Indices, cpt_flow = return_flow_indices(Aq_well_connect, Nt, Naq)
    #Indices = Indices - 1 #Fix the array starts at zero problem
    
    # List of boreholes that intersect each well
    WellsInterAq = return_boreholes_list(Aq_well_connect, Naq)
    
    # Computation of the integrals required in the linear system
    IntijI = integral_computation(Aq_well_connect,vect_t,L,Trans_aq,
        Stora_aq,R,Naq)

    # Definition of the matrices and the second member of the system
    Aglobal = np.zeros([cpt_flow+1, cpt_flow+1])
    B = np.zeros([cpt_flow, 1])

    # Loop over wells
    # pdb.set_trace()   
    
    for i in range(len(Aq_well_connect)):
        # Loop over the connected aquifers
        for I in range(len(Aq_well_connect[i])):
            # Definition of the aquifer numbering for the studied borehole
            num_aqi = Aq_well_connect[i][I]
            betaiI = betaI[i,I]
            # Loop over time
            for k in range(len(vect_t)):
                # Matrix value related to the studied borehole aquifer
                indexI1 = Indices[i,num_aqi,k] # well i and aquifer I
                #pdb.set_trace()
                Aglobal[indexI1,indexI1] = 1
                # Flow at the top of the pumped borehole
                if i==0 and I==0:
                    # For the pumping well: flow = pumped flow
                    if k*delta_t <= pump_time:
                        B[indexI1] = q_pumped
                    # General case
                else:
                    cpt_well = 0
                    while True:
#                   for cpt_well in range(len(WellsInterAq[num_aqi])): 
                       j = WellsInterAq[num_aqi][cpt_well]
                       # Number of wells intersected by borehole j
                       nb_aqj = len(Aq_well_connect[j])
                       # Convolution product in time
                       #for l in range(k):
                       l = 0    
                       while True:
                           # Relationship with aquifer I
                           indexI2 = Indices[j,num_aqi,l]
                           # intI1ij: Gamma_I, k^ij (aquifer I)
                           intIij_value = IntijI[i,j,num_aqi,k-l]
                           # flow qI
                           if I==0 and i==1:
                               Aglobal[indexI1,indexI2] += intIij_value/delta_t
                               if l > 0:
                                   # Indices for the well j and aquifer I
                                   indexI3 = Indices[j,num_aqi,l-1]
                                   Aglobal[indexI1,indexI3] -= intIij_value/delta_t
                           else:
                               Aglobal[indexI1,indexI2] += 0.5*betaiI*intIij_value
                               if l > 0:
                                   indexI4 = Indices[j,num_aqi,l-1]
                                   Aglobal[indexI1,indexI4] += 0.5*betaiI*intIij_value
                           # flow qI + 1
                           # (matlab) index=find(Aq_well_connect(j).aq_connect==num_aqi,1);
                           index = indices(Aq_well_connect[j], lambda x: x==num_aqi)[0]
                           if index < nb_aqj:
                               num_aq_nextj = Aq_well_connect[j][index+1]
                               indexI2 = Indices[j,num_aq_nextj,l]
                               if I==0 and i==1:
                                   Aglobal[indexI1,indexI2] -= intIij_value/delta_t
                                   if l > 0:
                                       indexI3 = Indices[j,num_aq_nextj,l-1]
                                       Aglobal[indexI1,indexI3] += intIij_value/delta_t
                               else:
                                   Aglobal[indexI1,indexI2] -= 0.5*betaiI*intIij_value
                                   if l > 0:
                                       indexI4 = Indices[j,num_aq_nextj,l-1]
                                       Aglobal[indexI1,indexI4] -= 0.5*betaiI*intIij_value
                           l += 1
                           if (l >= k ): break
                       cpt_well += 1
                       if (cpt_well >= len(WellsInterAq[num_aqi])): break  
                # Connection to the boreholes intersecting the previous aquifer I'
                # When it is not a flow at the top of the borehole
                    if I > 0:
                        num_aq_previousi = Aq_well_connect[i][I-1]
                        # This loop needs to execute at
                        # for cpt_well in range(len(WellsInterAq[num_aq_previousi])):
                        cpt_well = 0
                        while True:
                            j = WellsInterAq[num_aq_previousi][cpt_well]
                            nb_aqj = len(Aq_well_connect[j])
                            # Convolution in time
                            # for l in range(k):
                            l = 0
                            while True:
                                # intI1ij: Gamma_I, k^ij (aquifer I-1)
                                intIij_value = IntijI[i,j,num_aq_previousi,k-l]
                                # flow qI - 1
                                indexI2 = Indices[j,num_aq_previousi,k-l]
                                Aglobal[indexI1,indexI2] -= 0.5*betaiI*intIij_value
                                if l > 0:
                                    indexI4 = Indices[j,num_aq_previousi,l-1]
                                    Aglobal[indexI1,indexI4] -= 0.5*betaiI*intIij_value
                                    # flow qI
                                    index_previous = indices(Aq_well_connect[j],
                                        lambda x: x==num_aq_previousi)[0]
                                    if index_previous < nb_aqj:
                                        num_aq_nextj = Aq_well_connect[j][index_previous]
                                        indexI2 = Indices[j,num_aq_nextj,l]
                                        Aglobal[indexI1,indexI2] += 0.5*betaiI*intIij_value
                                    if l > 0:
                                        indexI4 = Indices[j,num_aq_nextj,l-1]
                                        Aglobal[indexI1,indexI4] += 0.5*betaiI*intIij_value
                                l += 1
                                if (l > k): break
                            cpt_well += 1
                            if (cpt_well >= len(WellsInterAq[num_aq_previousi])): break                
    # Linear system solution
    pdb.set_trace()                            
    q = np.linalg.solve(Aglobal,B)
    # Storage
    # Flowvelocities in the pumped borehole
    Naq_pump = len(Aq_well_connect[0])
    q_pump = np.zeros([Naq_pump, Nt])
    for I in range(Naq_pump):
        q_pump[I,:] = q[I*Nt:I*Nt+Nt]
    # Flow velocities in the observation borehole
    q_obs = np.zeros([Naq, Nt])
    for i in range(Naq):
        q_obs[I,:] = q[Nt*Naq_pump+I*Nt:Nt*Naq_pump+I*Nt+Nt]
    return (q_pump, q_obs) 
    
# ------------------------------------------------------------------------------
    
def fun_G1(r,theta,t1,t2,xij,alpha):
    beta = ((xij-r*cos(theta))**2 + (r*sin(theta))**2)/(4*alpha)
    if t1==0:
        res = r*exp1(beta/t2)
    else:
#        res = r*(expn(1,beta/t2) - expn(1,beta/t1))    
        res = r*(exp1(beta/t2) - exp1(beta/t1))   
    return res

# ------------------------------------------------------------------------------

def fun_G2(r,t,alpha):
    c = 4 * alpha * t
    if t==0:
        res = 0.0
    else:
#        res = c*pi*(1 - exp(-r**2/c)) + pi*r**2*expn(1,r**2/c)
        res = c*pi*(1 - exp(-r**2/c)) + pi*r**2*exp1(r**2/c)
    return res

# ------------------------------------------------------------------------------
def integral_computation(Aq_well_connect,vect_t,L,Trans_aq,Stora_aq,R,Naq):
    ''' Function to compute the integrals required in the linear system'''
    
    # Variables
    nwell = len(Aq_well_connect)
    alpha = Trans_aq/Stora_aq # Ratio be transmisivity and storativity
    delta_t = vect_t[1] - vect_t[0] # Timestep
    
    # Arrays to store the computed integrals for the wells i,j and 
    # aquifer I (must pass shape to zeros as a tuple)
    IntijI = np.zeros((nwell, nwell, Naq, vect_t.shape[0]))
    
    # J.E. NYQUIST: Note the original matlab code was written here to handle
    # more than two wells, but never tested for that.  I've simplified by
    # hardwiring the two well case.  This will need to be changed to add
    # additional wells.
    # For each well
    # pdb.set_trace()
    for i in range(nwell):
        for I in range(len(Aq_well_connect[i])):
            num_aq = Aq_well_connect[i][I]
            # Second loop on the wells to check if the aquifer intersects both
            for j in range(nwell):
                if (num_aq in Aq_well_connect[j]):
                    if i == j:
                        for k in range(vect_t.shape[0]): # Loop over time
                            IntijI[i,j,num_aq,k] = \
                                (fun_G2(R[j],vect_t[k],alpha[I])
                                - fun_G2(R[j],vect_t[k]-delta_t,alpha[I]))/ \
                                (4*pi*Trans_aq[I])
                            #pdb.set_trace()      
                    else:
                        xij = return_distance_wells(i,j,L)
                        for k in range(vect_t.shape[0]):
                            t1 = vect_t[k] - delta_t
                            t2 = vect_t[k]         
                            IntijI[i,j,num_aq,k] = dblquad(lambda r, theta: 
                                fun_G1(r,theta,t1,t2,xij,alpha[I]), 
                                0, 2*pi, lambda r: 0, lambda r: R[j])[0]
                            IntijI[i,j,num_aq,k] /= 4*pi*Trans_aq[I]
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
    #q_pump, q_obs = flow_velocity_computation(p)
    q_pump, q_obs = flow_velocity_computation(p.vect_t,p.Trans_aq,p.Stora_aq,
                    p.dist_aq,p.R,p.L,p.q_pumped,p.Aq_well_connect,p.pump_time,p.Naq)
    print(type(q_pump))
    pass
    pass
 
    #
    
    # 2.2 Conversion from [m/s] to [L/min]
    
    # 3. PLOT RESULTS FOR FLOW VELOCITY COMPARISON
