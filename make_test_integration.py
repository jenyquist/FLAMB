# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
from __future__ import division
from collections import namedtuple
from math import *
from scipy.special import expn
from scipy.integrate import dblquad
import numpy as np

        
def fun_G1(r,theta,t1,t2,xij,alpha):
    beta = ((xij-r*cos(theta))**2 + (r*sin(theta))**2)/(4*alpha)
    if t1==0:
        res = r*expn(1,beta/t2)
    else:
        res = r*(expn(1,beta/t2) - expn(1,beta/t1))      
    return res


# Test double quadrature on fun_G1
t1 = 0.0
t2 = 24.0
xij = -21.0
alpha = 2.0833
R = 0.0375

ans = dblquad(lambda r, theta: fun_G1(r,theta,t1,t2,xij,alpha), 0, 2*pi, lambda r: 0, lambda r: R)[0]

# quad2d(@(r_,theta_)fun_G1(r_,theta_,vect_t(k)-delta_t,vect_t(k),xij,alpha(I)),0,R(j),0,2*pi)/(4*pi*Trans_aq(I));

# Works!