# Python Finite-Difference-Code

The Python Finite-Difference code is tested with Python 2.7 and 3.4. The modules *numpy* and and *matplotlib* are required.

There are to versions of each file. The one without `_fast` in the file name is identical to the Matlab version. The `_fast` versions use a vectorization of the loops over the grid points, which will speed up the code significantly.
In Matlab this vectorization is done automatically in the background, therefore, no speed up could be archived by adding this to the Matlab code.  
