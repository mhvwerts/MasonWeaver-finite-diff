#!/bin/env python3
# -*- coding: utf-8 -*-

"""
Finite-difference solver for the Mason-Weaver equation.

The Mason-Weaver equation describes the sedimentation of Brownian particles,
i.e. it models sedimentation in the presence of diffusion. The Mason-Weaver
equation is a partial differential equation describing the time-evolution of
the particle concentration profile.

The present solver is based on the Crank-Nicolson finite-difference scheme
described in:

[1] Midelet, J.; El-Sagheer, A. H.; Brown, T.; Kanaras, A. G.;
    Werts, M. H. V. "The Sedimentation of Colloidal Nanoparticles in 
    Solution and Its Study Using Quantitative Digital Photography".
    Part. & Part. Syst. Charact. 2017, 34, 1700095. 
    https://doi.org/10.1002/ppsc.201700095 
    

by M. H. V. Werts, 2015-2020
Centre National de la Recherche Scientifique, France
Ecole normale supérieure de Rennes, lab. SATIE


Published under the CeCILL v2.1 license
    https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt


NO WARRANTIES. USE AT YOUR OWN RISK AND FOR THE GOOD OF MANKIND.
"""



import numpy as np
from scipy import sparse
from scipy.sparse import linalg



class aMWsoln:
    """Solution to the adimensional Mason-Weaver equation.
    
    The Mason-Weaver equation is solved with the finite-difference scheme 
    published in Ref. [1]. The entire finite-difference solution is calculated
    using the supplied parameters, immediately upon instantiation of the class.
    The finite-difference solution is then stored in a 2D Numpy array for
    further analysis.
    
    The solution is calculated on an exponentially expanding time grid, and an
    evenly spaced space grid.    
    
    Args:
        zeta_max (float): adimensional height of the cell
        tau_end (float): calculate solution from tau=0 to tau=tau_end
        J (int): number of space points minus 1 on (linear) grid
        N (int): number of time points minus 1 on (exponential) grid
        
        
    References:
        
    [1] Midelet, J.; El-Sagheer, A. H.; Brown, T.; Kanaras, A. G.;
        Werts, M. H. V. "The Sedimentation of Colloidal Nanoparticles in 
        Solution and Its Study Using Quantitative Digital Photography".
        Part. & Part. Syst. Charact. 2017, 34, 1700095. 
        https://doi.org/10.1002/ppsc.201700095 

    """
    def __init__(self, zeta_max = 10., tau_end = 10.,
                 J = 500, N = 100):
        # set calculation parameters
        self.zeta_max = float(zeta_max) # "height" of cell
        self.tau_end = float(tau_end) # end tau
        self.J = int(J) # number of space points
        self.N = int(N) # number of time points
        c_init = np.ones(self.J+1) # IC
        # define time and space grids
        self.k_tau = np.log(self.tau_end+1)/N
        self.tau = -1.0 + np.exp(self.k_tau*np.arange(0,self.N+1))
        self.dltzeta = self.zeta_max/self.J
        self.zeta=np.linspace(0,self.zeta_max,self.J+1,endpoint=True)
        # create c array for storing result
        self.c = np.zeros((self.J+1,self.N+1))
        self.c[:,0] = c_init
        # prepare
        c_n = c_init
        dltzeta=self.dltzeta
        tau=self.tau
        # loop generating c^{n+1} from c^{n} starting from c^0
        for n in range(0,N):
            # (re-)calculate alpha, beta, gamma
            # recalculation is actually only necessary for gamma
            alpha = 1/(2*(dltzeta**2))
            beta = 1/(4*dltzeta)
            gamma = 1/(tau[n+1]-tau[n])
            
            # construct LEFT HAND SIDE of matrix equation
            # finite difference diagonals
            ldiagelem = -alpha + beta
            cdiagelem = gamma + 2*alpha
            rdiagelem = -alpha - beta
            # boundary conditions LHS
            cstart = -2*beta + 0.25
            rstart = 2*beta + 0.25
            lend = -2*beta + 0.25
            cend = 2*beta + 0.25
            # create LHS tridiagonal matrix
            LHSmat = self._tridiamatrix(ldiagelem,cdiagelem,rdiagelem,
                                  cstart,rstart,lend,cend)
            
            # construct RIGHT HAND SIDE
            # finite difference diagonals
            ldiagelem = alpha - beta
            cdiagelem = gamma - 2*alpha
            rdiagelem = alpha + beta
            # boundary conditions RHS
            cstart = 2*beta - 0.25
            rstart = -2*beta - 0.25
            lend = 2*beta - 0.25
            cend = -2*beta - 0.25
            # create RLHS tridiagonal matrix
            RHSmat = self._tridiamatrix(ldiagelem,cdiagelem,rdiagelem,
                                  cstart,rstart,lend,cend)
            # construct RHS vector
            RHSvec = RHSmat * c_n # c contains concentration profile
            
            # SOLVE the matrix equation, giving c^{n+1} (denoted c_next)
            c_next = linalg.spsolve(LHSmat,RHSvec)
            # store c_{n+1} into solution matrix
            self.c[:,n+1] = c_next 
            # We swap the vectors containing c_n and c_next
            # in the next cycle c_next is then replaced based
            # on the new c_n.
            # Using the intermediate cswap reference is necessary, in order
            # to avoid that c_n and c_next end up referring to the same object
            cswap = c_n
            c_n = c_next
            c_next = cswap

    def _tridiamatrix(self,ldiagelem,cdiagelem,rdiagelem,
                 cstart,rstart,lend,cend):
        """utility function that generates a sparse tridiagonal 
        matrix from given elements"""
        ldiag = np.empty(self.J+1)
        cdiag = np.empty(self.J+1)
        rdiag = np.empty(self.J+1)
        ldiag.fill(ldiagelem)
        cdiag.fill(cdiagelem)
        rdiag.fill(rdiagelem)
        cdiag[0]=cstart
        rdiag[1]=rstart
        ldiag[-2]=lend
        cdiag[-1]=cend
        diag=[ldiag,cdiag,rdiag]
        N = self.J+1 # N=len(self.cdiag)
        return sparse.spdiags(diag,[-1,0,1],N,N,format="csr")
        
    def c_prof(self, tau_x):
        """Extract concentration profile at tau=tau_x (with interpolation)

        The interpolation will be done linearly, since we expect the time
        grid (more than) sufficiently dense        
        """
        # convenient references
        tau = self.tau
        c = self.c
        # find index (the point will lie between idx and idx-1)
        idx = np.searchsorted(tau, tau_x, side="left")
        if ((idx<0) or idx>=len(tau)):
            raise Exception('tau_x out of tau range')
        frac=(tau_x - tau[idx-1])/(tau[idx] - tau[idx-1])
        c_interp = (1.-frac)*c[:,idx-1] + frac*c[:,idx]
        return c_interp



class MWsoln(aMWsoln):
    """Solution to the Mason-Weaver equation, using real-world (SI) units.
    
    The Mason-Weaver equation is solved with the finite-difference scheme 
    published in Ref. [1]. The entire finite-difference solution is calculated
    using the supplied parameters, immediately upon instantiation of the class.
    The finite-difference solution is then stored in a 2D Numpy array for
    further analysis.
    
    The solution is calculated on an exponentially expanding time grid, and an
    evenly spaced space grid. There are J+1 space points and N+1 time points.
    The corresponding positions and times are in the vectors z and t, 
    respectively, that are properties of the object.    

    Args:
        z_max (float): height of the cell  [m]
        t_end (float): calculate solution from t=0 to t=t_end  [s]
        J (int): number of space points minus 1 on (linear) grid
        N (int): number of time points minus 1 on (exponentially expanding)
                 grid
        D (float): Fickian diffusion coefficient [m^2 s^{-1}]
        sg (float): product of sedimentation coefficient and gravitational
                    acceleration, overall: [m s^{-1}]        
        
    References:
        
    [1] Midelet, J.; El-Sagheer, A. H.; Brown, T.; Kanaras, A. G.;
        Werts, M. H. V. "The Sedimentation of Colloidal Nanoparticles in 
        Solution and Its Study Using Quantitative Digital Photography".
        Part. & Part. Syst. Charact. 2017, 34, 1700095. 
        https://doi.org/10.1002/ppsc.201700095 

    """
    def __init__(self, z_max, t_end, J, N, D, sg):
        self.z_max = float(z_max)
        self.t_end = float(t_end)
        self.D = float(D)
        self.sg = float(sg)
        self.z0 = self.D/self.sg
        self.t0 = self.D/(self.sg**2)
        zeta_max = self.z_max/self.z0
        tau_end = self.t_end/self.t0
        super(MWsoln, self).__init__(zeta_max, tau_end, J, N)
        self.t = self.tau * self.t0
        self.z = self.zeta * self.z0
        
    def c_prof(self, t_x):
        """Extract concentration profile at t=t_x (with interpolation)

        If t_x is between two calculated time points, a linear interpolation 
        between the two closest calculated time points will be done, expecting
        the time grid to be more than sufficiently dense.     
       

        Parameters
        ----------
        t_x : float
            Time t=t_x for which Mason-Weaver concentration profile is requested.

        Returns
        -------
        np.array 
            concentration profile at time t=t_x.

        """
        tau_x = t_x/self.t0
        return super(MWsoln, self).c_prof(tau_x)
    
    

# Example code
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    Npts = 500

    soln = MWsoln(z_max=0.01, t_end=3.6E6,
                  J=1000, N=Npts,
                  D=4.32E-12, 
                  sg=2.47E-9*9.81) # for 60 nm gold nanoparticle at 4C
    
    plt.figure(1)
    plt.clf()
    plt.title('Sedimentation 60nm diam. gold NP in water at 4°C')
    for i in range(0, Npts, 50):
        plt.plot(soln.z[:-1]*100, soln.c[:-1,i],
                 label='t={:.0f} h '.format(soln.t[i]/3600))
    plt.ylim(ymin = 0.0, ymax = 3.0)
    plt.xlabel('height / cm')
    plt.ylabel('rel. concentration')
    plt.legend()
    
    # time behaviour
    plt.figure(2)
    plt.clf()
    plt.title('Concentration at bottom of cell')
    plt.plot(soln.t/3600, soln.c[1,:].T)
    plt.ylabel('rel. concentration')
    plt.xlabel('time / h')
    
    # convergence to analytical steady-state
    plt.figure(3)
    plt.clf()
    plt.title('Sedimentation equilibrium gradient')
    B = soln.z_max / (soln.z0*(1-np.exp(-soln.z_max/soln.z0)))
    plt.plot(soln.z*100, soln.c[:,-1],
             label = 'finite diff.')
    plt.plot(soln.z*100, B * np.exp(-soln.z/soln.z0),
             label = 'analytic')
    plt.xlabel('height / cm')
    plt.ylabel('rel. concentration')
    plt.legend()
    
    plt.show()
    
    
