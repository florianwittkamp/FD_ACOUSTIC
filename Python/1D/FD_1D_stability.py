## FD_1D_stability.py 1-D FD stability calculation
# GNU General Public License v3.0
#
# Author: Florian Wittkamp 2016
#
# Calculate stability limit of FD-simulations
#
# Stability limit is calculated in terms of the CFL-number, which is
# defined as: CFL=v_(max)*DT/DX
# You get the maximum DT by DT=CFL*DX/v_(max)
#
# Theory:
# Fei, X., & Xiaohong, T. (2006).
# Stability and numerical dispersion analysis of a fourth-order accurate FDTD method. Antennas and Propagation,
# IEEE Transactions on, 54(9), 2525-2530.

## Initialisation
import numpy as np
from FD_taylor_coeff_func import coeff
print(" ")
print("Starting FD_1D_stability")

## Input Parameter
spatial_order=4

## Calculate spatial derivation impact
# Assumption: Spatial sampling is at the Nyquist condition
sum_fd_stencil=np.sum(np.abs(coeff(spatial_order)))
theta=np.arange(0,2*np.pi,0.01)
print("You choose a spatial order of ",spatial_order)

## Temporal 2 order (Leapfrog)
f=-1j*(((np.exp(-1j*theta))**(2.5)-(np.exp(-1j*theta))**1.5))
f2=(2*sum_fd_stencil)
f_ges=f/f2
print("2. Temporal Order has the stability limit CFL=",np.max(f_ges.real))

## Temporal 3 order (Adams-Bashforth method)
f=-1j*(((np.exp(-1j*theta))**2.5-(np.exp(-1j*theta))**1.5)*1);
f2=(2*sum_fd_stencil*(25.0/24.0*(np.exp(-1j*theta))**2
    -1./12.*(np.exp(-1j*theta))**1+1./24.))
f_ges=f/f2
print("3. Temporal Order (ABS-method) has the stability limit CFL=",np.max(f_ges.real))

## Temporal 4 order (Adams-Bashforth method)
f=-1j*(((np.exp(-1j*theta))**3.5-(np.exp(-1j*theta))**2.5))
f2=(2*sum_fd_stencil*(13./12.*(np.exp(-1j*theta))**3
    -5./24.*(np.exp(-1j*theta))**2+1./6.*(np.exp(-1j*theta))**1-1./24.))
f_ges=f/f2
print("4. Temporal Order (ABS-method) has the stability limit CFL=",np.max(f_ges.real))

print(" ")
