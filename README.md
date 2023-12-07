# Magnetorotational Instability （MRI) 
The magnetorotational instability (MRI) is a very important instability in astrophysics. It is a fluid instability that causes an accretion disk orbiting a massive central object to become trubulent.  
  
The MRI happens when a group of plasma rotating around a central object with big mass (like black hole) and will destroy the structure of the plasma.We are going to use Python 3 to simulate the MRI phenomenon in two dimension. We mainly refer to the paper written by Riquelme[1].  
  
We use Particle-In-Cell (PIC) method which divide the box into many grids(cells) and put several particles in each grid. We will only record the position of the particles. For other variables like velocity, currency, electric field, and magnetic field, we assume that they are not attached the particles but the grids. 

# Equations for the MRI Tested

Our initial equations for the MRI in a shearing coordinate frame are (Riquelme, 2012)
$$\frac{\partial\vec{B}}{\partial t} = -\nabla \times \vec{E} - sB_x\hat{y}$$
$$\frac{\partial\vec{E}}{\partial t} = \nabla \times \vec{B} -4\pi\vec{J}- sE_x\hat{y}$$
$$\frac{d\vec{p}}{dt} = 2\omega_0p_y\hat{x} - \frac{1}{2}\omega_0p_x\hat{y} + q(\vec{E}+\frac{\vec{u}}{c}\times\vec{B})$$

We normalized the electric and magnetic fields to initial background z component of the magnetic field, the time to the orbital frequency, and the velocity to the Alfvén speed. 

# Two Approaches

We have two simulations to explore the MRI.

Code 1. SIVP Code: This code uses the built-in SciPy.Integrate routine solve_ivp to directly integrate the equations for the magnetorotational instability.
More details about inputs and code options are described in the file SIVPSimulation.ipynb, and some preliminary results are described in SIVPSimulation.ipynb.
This code produces results in SIVPAnalysis.ipynb similar to those found by Riquelme et al., but more careful study with improved grid resolution is needed to determine if our results are consistent with the present literature.
Current results are limited by instability of a vectorized scheme to do PIC weighting.

Code 2. Central Time Code:
This code imports and uses scipy, numpy and math.
The function is initialized by: 
  PIC = Particle_in_cell_2D(Lx=0.1,Lz=1,Nx=10,Nz=100,t_max=10, dt=0.01, dt_out = 0.1, alpha=.99999,D=1)
  
Its inputs are:
  Lx - Length of grid in x
  Lz - Length of grid in z
  Nx - Number of cells in x
  Nz - Number of cells in z
  t_max - Maximum time to run
  dt_out - output times, rate at which to save fields and particle positions
  alpha - effectively a momentum and rejection term. The intent was to provide small smoothing of the field updates.
  D - a diffusive factor to the laplace smoothing

  Lx, Lz, Nx, Nz, t_max are mandatory inputs. dt_out defaults to 0.1, alpha and D default to 1. 
  
Its used as:  
  x, z, pos_x, pos_z, Ux, Uy, Uz, Bx, By, Bz, Ex, Ey, Ez = PIC.solve()

The outputs are as follows:
  x - 1D array of x dimension
  z - 1D array of z dimension
  pos_x - 2D list of particle positions in x for each times stamp of dt_out + initial condition
        - time is the first axis; Ex: posx[0] returns array of x particle positions at initial condition
  poz_z - ^ but for z
  Ux, Uy, Uz - 3D list of particle velocities at each dt_out + IC. Time is first axis, x is second, z is third.
    Example: Ux[0][0,:] returns the z dependant Ux at t = 0, x = 0. 
  Bx, By, Bz - ^ but for magnetic field
  Ex, Ey, Ez - ^ but for electric fields

*** As is, the code is not stable. ***
Without smoothing, i.e. Alpha = 1 and D=1, the code can run to about t = 1.7, 
before we see massive particle runaway due to sharp gradients in the fields and velocities. 
For a small amount of smoothing, Alpha = 0.99999 and D = 1, we reach longer times of about t = 10 
before we observe the same breakdown. I think this isn't necessarily a mark of increased stability,
but a slowdown of the dynamics of the problem. 

We wanted to prevent particles from traveling across more than one cell in one time step by implementing
adaptive time scaling. It works really well, but unfortunately when the code breaks down we see velocities
that are just too sharp and the time steps drop drastically. Additionally, we could not decrease grid size
because immediately out time step would plummet since the velocities were relatively large and we couldn't
contain one particle within one grid cell. 

All this to say, the input values provided above allow for simulation to times just before the breakdown 
of out system. We tried very much to find the cause of this but were unsuccessful, none of our fixes worked. 

It is possible that normalizations of the fields weren't done correctly despite rigorous checks. If this work
was to be continued, I think it would be good to check the size of each term. Analysis has shown that without
the current density terms, that is curl B - J = 0, led to stable results. Then it must be in the current density
term / particle position terms which have caused the problems. 

Plots of the simulation for alpha = 1 and alpha = 0.99999 respectively. In Alpha 1 we show the moments before
breakdown of the simulation with Ux. 

# Example and Results

Alpha = 1
![image](https://github.com/ehansen99/Python-Kinetic-Magnetorotational-Instability/assets/143833778/4cdf8a66-7d87-48a7-bed0-cbf0ca7cc257)

![image](https://github.com/ehansen99/Python-Kinetic-Magnetorotational-Instability/assets/143833778/3ef7a65c-007b-4093-9950-5ce6c9290bc0)

Alpha = 0.99999
![image1](https://github.com/ehansen99/Python-Kinetic-Magnetorotational-Instability/assets/107236110/4b2c269b-906f-43c2-9cbf-11cc0a50266d)

# Requirements
Packages we use:  
&bull; numpy  
&bull; matplotlib  
&bull; scipy  
&bull; jupyter

The SIVP simulation specifically makes use of
* scipy.sparse
* scipy.integrate
  
and further uses the external package
* opt_einsum
  
for optimum weighting array contractions.


# References

S. A. Balbus and J. F. Hawley. A Powerful Local Shear Instability in Weakly
Magnetized Disks. I. Linear Analysis., 376:214, July 1991.

C. K. Birdsall and A. B. Langdon. Plasma physics via computer simulation /
.K. Birdsall, A.B. Langdon. Series in Plasma Physics. Adam Hilger, Bristol,
1991.

M. A. Riquelme, E. Quataert, P. Sharma, and A. Spitkovsky.
Local two-dimensional particle-in-cell simulations of the collisionless magne-
torotational instability. The Astrophysical journal, 755(1):1–20, 2012.

Daniel G. A. Smith and Johnnie Gray, opt_einsum - A Python package for optimizing contraction order for einsum-like expressions. Journal of Open Source Software, 2018, 3(26), 753.
