# ------------------------------------------------------------------------------
from __future__ import division
from collections import namedtuple
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
    
    # Check parameters
    assert len(dist_aq0) == Naq, 'Number of aquifers wrong in pumping well.'
    assert len(dist_aq1) == Naq, 'Number of aquifers wrong in observation well'
    assert len(Stora_aq) == Naq, 'Number of storativities is wrong'
    
    params = namedtuple('params',['L', 'R', 'pump_time', 'Trans_aq', 'Stora_aq',
                    'vect_t', 'dist_aq', 'q_pumped', 'Naq', 'Aq_well_connect'])
    p = params(L,R,pump_time,Trans_aq,Stora_aq,vect_t,dist_aq,q_pumped,Naq,Aq_well_connect)
    return p

    # return L,R,pump_time,Trans_aq,Stora_aq,vect_t,dist_aq,q_pumped,Naq,Aq_well_connect

# ------------------------------------------------------------------------------
    
def return_flow_indices(Aq_well_connect, Nt, Naq):
    ''' Return the number for the flow in the linear system for a given 
    well i and aquifer num_aq'''
        
    Indices = np.zeros((len(Aq_well_connect), Naq, Nt))
    cpt_flow = 0  # Python arrays start at zero
    
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

    def IntijI = IntergralComputation(p):
        ''' Function to compute the integrals required in the linear system'''
        
        # Variables
        nwell = len(p.Aq_well_connect)
        alpha = p.Trans_aq/p.Stora_aq
        delta_t = p.vect_t[1] - p.vect_t[0]
        
        # Arrays to store the computed integrals for the wells i,j and 
        # aquifer I
        IntijI = np.zeros(nwell, nwell, p.Naq, len(p.vect_t))
        
        # For each well
        for i in range(nwell):
            # Loop over the aquifers connected to this well
            for I in range(len(p.Aq_well_connect[i]):
                # Second loop on the wells to see if the well interests the
                # aquifer in question
                for j in range(nwell):
                    pass
        return IntijI
                    
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
    betaI = np.zeros(p.dist_aq.shape)
    betaI[0,:] = rho*g*p.R[0]**2/(8*mu*p.dist_aq[0,:]) # pumped well
    betaI[1,:] = rho*g*p.R[1]**2/(8*mu*p.dist_aq[1,:]) # observation well
    
    # 1.0 LINEAR SYSTEM DEFINITION
    
    # Flow indices
    Indices, cpt_flow = return_flow_indices(p.Aq_well_connect, Nt, p.Naq)
    
    # List of boreholes that intersect each well
    WellsInterAq = return_boreholes_list(p.Aq_well_connect, p.Naq)
    
    # Computation of the integrals required in the linear system
    IntijI = IntergralComputation(p)
    
    return (WellsInterAq)
    
  
    
    
    
    
    # ------------------------------------------------------------------------------