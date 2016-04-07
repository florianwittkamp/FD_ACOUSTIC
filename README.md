# Finite-Difference seismic wave simulation

This is a collection of scripts to simulate seismic wave propagation in 1-D and 2-D.
The wave propagation is based on the first-order acoustic wave equation in stress-velocity formulation (Virieux (1986)), which is solved by Finite-Differences. 

# Content

At the moment only 1-D versions of the FD wave simulation codes are available. The 2-D versions will be available in the future. 

The 1-D versions are available in the folder Matlab/1D for different spatial and temporal orders. 

Higher spatial orders are achived by a classical Taylor expansion.

For higher temporal orders there are two methods available:

1. Lax–Wendroff method (Theory: Dablain (1986))

2. Adams-Bashforth method (Theory: Bohlen & Wittkamp (2016))

To compare the effect of different orders of accuracy you can run the script FD_compare.m.

# Literatur

* Bohlen, T., & Wittkamp, F. (2016). Three-dimensional viscoelastic time-domain finite-difference seismic modelling using the staggered Adams–Bashforth time integrator. Geophysical Journal International, 204(3), 1781-1788.

* Dablain, M. A. (1986). The application of high-order differencing to the scalar wave equation. Geophysics, 51(1), 54-66.

* Virieux, J. (1986). P-SV wave propagation in heterogeneous media: Velocity-stress finite-difference method. Geophysics, 51(4), 889-901.


# Licence 
This collection is available under the GNU General Public License v3.0. See the LICENCE file for more information.
