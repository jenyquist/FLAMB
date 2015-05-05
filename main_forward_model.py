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
import xflow

# 1. PARAMETER DEFINITION
# All of the model parameters are set in this function
# p is a named tuple with the following parameters:
# L,R,pump_time,Trans_aq,Stora_aq,vect_t,dist_aq,q_pumped,Naq,Aq_well_connect
# Each can be accessed with the dot notation (eg. p.pump_time)
p = xflow.parameters()

# 2. FLOW VELOCITY DETERMINATION
# 2.1 Flow velocity computation from the developed forward model
# q_pump, q_obs = flow_velocity_computation(p)
WellsInterAq = xflow.flow_velocity_computation(p)
print(WellInterAq)
#

# 2.2 Conversion from [m/s] to [L/min]

# 3. PLOT RESULTS FOR FLOW VELOCITY COMPARISON