#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 16:31:22 2020

@author: arran

example:
    m*x''(t) + b*x'(t) + k*x(t)+a*(x(t))^3 = -m*g
    
    x' = y
    y' = -b/m*y - k/m*x - a/m*x**3 - g
    np.array([y, -(b*y + (k + a*x*x)*x) / m - g])

"""




# EQUATIONS 1 - 9 in Jenkins '11
# I think used the reduced ones, 14 - 19
# d/dx (D*U)  = e_dot + m_dot               (1)

# d/dx (D*U*U) =  D*((rho_a - rho)/ rho_0)*g*sin(alpha) - C_d*U*U

# d/dx (D*U*T ) = e_dot*T_a + m_dot*T_b - C**(1/2)*U*gamma_T*(T - T_b)

# d/dx (D*U*S) =  e_dot*S_a + m_dot*S_b - C**(1/2)*U*gamma_S*(S - S_b)

# rho = rho_0( 1 +beta_S*(S - S_0) - beta_T*(T - T_0))
# e_dot = E_0*U*np.sin(alpha)


# A = DU
# B = DUU
# C = DUT
# D = DUS

# A = D*U
# B = A*U
# C = A*T
# E = A*S 




# def func(y):
#     A,B,C,E,D,U,T,S = y.split(8)
    
#     dAdx  =   = e_dot + m_dot
#     dBdx = D*(rho_a - rho)/ rho_0*g*np.sin(alpha) - C_d*U*U
#     dCdx = e_dot*T_a + m_dot*T_b - C**(1/2)*U*gamma_T*(T - T_b)
#     dEdx = e_dot*S_a + m_dot*S_b - C**(1/2)*U*gamma_S*(S - S_b)
    
#     0 = rho - rho_0( 1 +beta_S*(S - S_0) - beta_T*(T - T_0))
#     0 = e_dot - E_0*U*np.sin(alpha)
    
#     0 = D*U - A
#     0 = A*U - B
#     0 = A*T - C
#     0 = A*S - E
    
    
#     return np.hstack(dAdx, dBdx, dCdx, dEdx)
    