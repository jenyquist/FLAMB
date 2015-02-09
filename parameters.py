from __future__ import division
from collections import namedtuple
import numpy as np
def parameters():
    '''This function is the place to set all of the parameters.'''
    
    #PARAMETER DEFINITIONS
    
    # 1.1 Borehole parameters
    
    # Distance between boreholes [meters]
    # (approximated as similar in each aquifer)
    L = 21
    
    # Radius of the pumped (i=0) and observation (i=1) boreholes [meters]
    R = [3.75e-2, 3.75e-3]
    
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
    Aq_well_connect = {'pump_well_connect':[2], 
                       'obs_well_connect': [1, 2, 3, 4, 5]}
                       
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