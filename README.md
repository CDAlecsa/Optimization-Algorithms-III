# Optimization-Algorithms-III
- optimizers for 1D and 2D objective functions (including Rosenbrock function). A comparison based upon the energy dissipation, phase portrait, decrease in the objective function
- this repository is based on the paper 'Asymptotic behaviour of a Strang splitting scheme for damped second-order gradient systems' 
by Cristian Daniel Alecsa, Titus Pinta & Adrian Viorel

- the optimization algorithms are presented in the attached '.pdf' file
- one cand modify the input parameters in the '.txt' file (in both folders)

Optimization Algorithms :
- Crank Nicolson (implicit method)
- Heun Method (explicit method)
- Strang Splitting (implicit method)
- Strang Splitting with Predictor-Corrector (explicit method)
- Polyak method (explicit method)
- Runge Kutta integrator

Objective functions :
- 1D
- linear 2D
- Rosenbrock
- a 2D sin type function that satisfies the PL inequality, from the article 'A geometric integration approach to smooth optimisation:
Foundations of the discrete gradient method' by M.J. Ehrhardt et. al.


The '.mw' files are made in Maple. These files contain the integrals from the RHS of the Strang splitting and Strang with Predictor-Corrector (for the 1D objective function) and the exact solution for the linear system of ODE's (in the 2D folder). Also, for Rosenbrock and 2D sin-type function, the '.mw' files contains the computed integrals that are used in '.py' files for the Newton methods for solving implicit systems.
