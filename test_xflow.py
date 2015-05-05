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
def test_parameters():
    # Check parameters
    p = xflow.parameters()
    assert len(p.dist_aq[0]) == p.Naq, 'Number of aquifers wrong in pumping well.'
    assert len(p.dist_aq[1]) == p.Naq, 'Number of aquifers wrong in observation well'
    assert len(p.Stora_aq) == p.Naq, 'Number of storativities is wrong'

def test_return_flow_indices():
    import pickle
    f = open("p.pkl","rb")
    p = pickle.load(f)
    Nt = len(p['vect_t'])
    (Indices, cpt_flow) = xflow.return_flow_indices(p['Aq_well_connect'], Nt, p['Naq'])
    assert Indices.shape == (2, 5, 100)

def test_integrate():
    import pickle
    f = open("p.pkl","rb")
    p = pickle.load(f)
    assert p['Naq'] == 5

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

def test_fun_G2():
    alpha = 2.0833
    t = 10.0
    r = 1.5
    answer_from_matlab = 28.614519953675380
    tol = 1e-8;
    assert abs(xflow.fun_G2(r,t,alpha) - answer_from_matlab) < tol

'''
Objects which have been dumped to a file with pickle.dump can be reread into a program by using the method pickle.load(file). pickle.load recognizes automatically, which format had been used for writing the data.
A simple example:

>>> cities = ["Paris", "Dijon","Lyon","Strasbourg"]
>>> fh = open("data.pkl","bw")
>>> pickle.dump(cities,fh)
>>> fh.close()

The file data.pkl can be read in again by Python in the same or another session or by a different program:

>>> import pickle
>>> f = open("data.pkl","rb")
>>> villes = pickle.load(f)
>>> print(villes)
['Paris', 'Dijon', 'Lyon', 'Strasbourg']
>>>
'''
