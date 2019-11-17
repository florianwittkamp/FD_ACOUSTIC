
# Finite-Difference seismic wave simulation

This is a collection of **Matlab** and **Python** scripts to simulate **seismic wave propagation** in 1-D and 2-D.
The wave propagation is based on the first-order acoustic wave equation in stress-velocity formulation (e.g. *Virieux (1986)*), which is solved by Finite-Differences on a staggered-grid.

### Content

In this repository 1-D and 2-D versions of Finite-Difference wave simulation codes are available in Matlab and Python.
The source code can be found in the directory `Matlab/`, `Python/` and `JupyterNotebook`, respectively.

Higher spatial-orders are achieved by a classical Taylor expansion.

For higher temporal-orders there are two methods available:

1. Lax–Wendroff method (only in 1-D). Theory: *Dablain (1986)*

2. Adams-Bashforth method. Theory: *Bohlen & Wittkamp (2016)*

To explore the influence of different orders of accuracy you can run the script `FD_1D_compare` or `FD_2D_compare`.

Moreover, in 1-D there scripts to calculate and plot the numerical dispersion as well as the numerical dissipation (Adams-Bashforth method) are provided. Up to now, this scripts are Matlab only. The underlying theory is given in *Bohlen & Wittkamp (2016)*.

### Literature

* *Bohlen, T., & Wittkamp, F.* (2016). Three-dimensional viscoelastic time-domain finite-difference seismic modeling using the staggered Adams–Bashforth time integrator. Geophysical Journal International, 204(3), 1781-1788.

* *Dablain, M. A.* (1986). The application of high-order differencing to the scalar wave equation. Geophysics, 51(1), 54-66.

* *Virieux, J.* (1986). P-SV wave propagation in heterogeneous media: Velocity-stress finite-difference method. Geophysics, 51(4), 889-901.


### Licence

This collection is available under the GNU General Public License v3.0. See the `LICENCE` file for more information.
