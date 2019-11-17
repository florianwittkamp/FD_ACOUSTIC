## FD_taylor_coeff_func.py calculate Taylor coefficient
# GNU General Public License v3.0
#
# Author: Florian Wittkamp 2016
#
# Calculate the Taylor coefficient in an arbitrary order for the
# second-order derivative.
#
# Usage:
# Lets say you want to calculate the 4th order accurate second-order
# FD-stencil. Then you have to set order=4 and the result would be used
# as follow:
# p_x = 1/DH * ( coeff(1) * (p(x+1)-p(x)) + coeff(2) * (p(x+2)-p(x-1)) )
# where p_x is the derivative.
import numpy as np
def coeff(order):
    ## Check some conditions
    if np.int(order)%2!=0:
        print("Error: coeff")
        print("Order has to be an integer multiple of 2!")
        return
    if order==2:
        print("Error: coeff")
        print("Order has to be at least 4!")
        return
    ## Calculation
    c=np.transpose(np.hstack((1, np.zeros(np.int(order/2)-1))))
    M=np.zeros((np.int(order/2),np.int(order/2)))
    # Condition 1: \sum^{N/2}_{k=1} b_k(2k-1)=1
    for n in range(1,np.int(order/2+1)):
        M[0,n-1]=(2*n-1)
    # Condition 2:  \sum^{N/2}_{k=1} b_k(2k-1)^(2j-1)=0; j=2,3...N/2
    for j in range(2,np.int(order/2+1)):
        for n in range(1,np.int(order/2+1)):
            M[j-1,n-1]=(2*n-1)**(2*j-1)
    coeff=np.transpose(np.dot(np.linalg.inv(M),c))
    return(coeff)
