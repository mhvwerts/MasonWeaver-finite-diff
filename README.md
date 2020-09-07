# MasonWeaver-finite-diff: Finite-difference solver for the Mason-Weaver equation.

The Mason-Weaver equation describes the sedimentation of Brownian particles,
*i.e.* it models sedimentation in the presence of diffusion. The Mason-Weaver
equation is a partial differential equation describing the time-evolution of
the particle concentration profile.

The present solver is based on the Crank-Nicolson finite-difference scheme
described in:

[1] Midelet, J.; El-Sagheer, A. H.; Brown, T.; Kanaras, A. G.;
    Werts, M. H. V. "The Sedimentation of Colloidal Nanoparticles in 
    Solution and Its Study Using Quantitative Digital Photography".
    Part. & Part. Syst. Charact. 2017, 34, 1700095. 
    [doi:10.1002/ppsc.201700095](https://doi.org/10.1002/ppsc.201700095) 


## Installation and usage

There is no specific installation. The program will run on Python 3 and depends on ``numpy`` and ``scipy``. It is actually a Python module from which the object class ``MWsoln`` can be imported. When creating an instance of this class, a finite-difference solution is immediately calculated (using the supplied parameters) and stored in an array that is accessible as a property of the class.

Inside the module, there is example code that demonstrates how to run ('instantiate') the solver and extract the results.


## To-do's

- Further documentation and examples
- Provide test code





