# Testing functions for xflow.py
# Since the matlab version works, many of these tests simply compare the output here with
# what matlab returns. 
# A sample set of parameters has been pickled on disk as p.pkl
# These parameters were pickled as an ordered dictionary as pickle has problems with named tuples.
# Tests can executed from a terminal shell with the command:
# py.test -s -v
#
# J. E. Nyquist, April 28, 2015

import xflow
import pdb
#pdb.set_trace()

def test_parameters():
    # Check parameters
    p = xflow.parameters()
    assert len(p.dist_aq[0]) == p.Naq, 'Number of aquifers wrong in pumping well.'
    assert len(p.dist_aq[1]) == p.Naq, 'Number of aquifers wrong in observation well'
    assert len(p.Stora_aq) == p.Naq, 'Number of storativities is wrong'

def test_return_boreholes_list():
    import pickle
    f = open("p.pkl","rb")
    p = pickle.load(f)
    WellsInterAq = xflow.return_boreholes_list(p['Aq_well_connect'], p['Naq'])
    assert WellsInterAq[1] == [0,1]

def test_return_conductivities():
    import pickle
    f = open("p.pkl","rb")
    p = pickle.load(f)
    betaI = xflow.return_conductivities(p['R'], p['dist_aq'])
    print(betaI)
    assert abs(betaI[0,4] - 66.3236) < 1e-4
    assert abs(betaI[1,4] - 74.9745) < 1e-4

def test_return_flow_indices():
    import pickle
    f = open("p.pkl","rb")
    p = pickle.load(f)
    Nt = len(p['vect_t'])
    (Indices, cpt_flow) = xflow.return_flow_indices(p['Aq_well_connect'], Nt, p['Naq'])
    assert Indices.shape == (2, 5, 100)
    #assert Indices[0,1,-1] == 100
    
    Aq_well_connect = [[1], [0, 1, 2, 3, 4]]
    Nt = 2
    Naq = 5
    Indices, cpt_flow = xflow.return_flow_indices(Aq_well_connect, Nt, Naq)
    assert cpt_flow == 11

def test_integrate():
    import pickle
    f = open("p.pkl","rb")
    p = pickle.load(f)
    assert p['Naq'] == 5

def test_fun_G1():
    alpha = 2.0833
    t1 = 10.0
    t2 = 12.5
    theta = 1.0
    r = 1.5
    xij = -21
    answer_from_matlab = 0.002082057384242
    tol = 1e-8;
    assert abs(xflow.fun_G1(r,theta,t1,t2,xij,alpha) - answer_from_matlab) < tol

def test_fun_G2():
    alpha = 2.0833
    t = 10.0
    r = 1.5
    answer_from_matlab = 28.614519953675380
    tol = 1e-8;
    assert abs(xflow.fun_G2(r,t,alpha) - answer_from_matlab) < tol
    alpha = 2.314814814814814
    t = 60.0
    r = 0.0375
    answer_from_matlab = 0.058799936084614
    tol = 1e-8;
    assert abs(xflow.fun_G2(r,t,alpha) - answer_from_matlab) < tol

def test_double_quadrature_on_fun_G1():
    from scipy.integrate import dblquad
    from math import pi
    t1 = 0.0
    t2 = 24.0
    xij = -21.0
    alpha = 2.0833
    R = 0.0375
    answer_from_matlab = 1.6319e-04
    ans = dblquad(lambda r, theta: xflow.fun_G1(r,theta,t1,t2,xij,alpha), 
            0, 2*pi, lambda r: 0, lambda r: R)[0]
    assert abs(ans - answer_from_matlab) < 1e-8

def test_integral_computation():
    import pickle
    import numpy as np
    import scipy.io
    f = open("p.pkl","rb")
    p = pickle.load(f)
    Aq_well_connect = p['Aq_well_connect']
    vect_t = p['vect_t']
    L = p['L']
    Trans_aq = p['Trans_aq']
    Stora_aq = p['Stora_aq']
    R = p['R']
    Naq = p['Naq']
    IntijI = xflow.integral_computation(Aq_well_connect,vect_t,L,Trans_aq,
            Stora_aq,R,Naq)
    assert abs(IntijI[0,1,1,4] - 2.2986) < 1e-4
    
    # A more complete test
    mat = scipy.io.loadmat('IntijI.mat')
    Aq_well_connect = [[1], [0, 1, 2, 3, 4]]
    vect_t = np.array([60, 120])
    L = 21
    R = np.array([0.0375, 0.0375])
    Trans_aq = np.array([0.208333333333333e-4, 0.231481481481481e-4, 
                         0.231481481481481e-4, 0.231481481481481e-5, 
                         0.231481481481481e-5,])
    Stora_aq = np.array([1e-5, 1e-5, 1e-5, 1e-5, 1e-5])
    IntijI = xflow.integral_computation(Aq_well_connect,vect_t,L,Trans_aq,
            Stora_aq,R,Naq)
    assert np.allclose(mat['IntijI'], IntijI)

def test_indices():
    a = [1, 2, 3, 1, 2, 3, 1, 2, 3]
    inds = xflow.indices(a, lambda x: x > 2)[0] #[0] being the 2nd matlab argument
    assert inds == 2

def test_xflow():
    import numpy as np
    tol = 1e-8;
    vect_t = np.array([60, 120])
    Trans_aq = np.array([0.208333333333333e-4, 0.231481481481481e-4, 
                         0.231481481481481e-4, 0.231481481481481e-4, 0.231481481481481e-4,])
    Stora_aq = np.array([1e-5, 1e-5, 1e-5, 1e-5, 1e-5])
    dist_aq = np.array([[26, 30, 40, 20, 26], [26, 26, 40, 24, 23]])
    R = np.array([0.0375, 0.0375])
    L = 21
    q_pumped = 0.067151596729736
    Aq_well_connect = [[1], [0, 1, 2, 3, 4]]
    pump_time = 120
    Naq = 5
    (q_pump, q_obs) = xflow.flow_velocity_computation(vect_t,Trans_aq,Stora_aq,dist_aq,R,L,q_pumped,Aq_well_connect,pump_time,Naq)
    qp = 0.0671515967297360
    qobs = 7.87931805456793e-05
    assert abs(qp - q_pump[1]) < tol
    assert abs(qobs - q_obs[4][1]) < tol
    