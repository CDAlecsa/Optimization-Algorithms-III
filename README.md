# Optimization-Algorithms-III
- optimizers for 1D and 2D objective functions (including Rosenbrock function). A comparison based upon the energy dissipation, phase portrait, decrease in the objective function
- this repository is based on the paper 'Asymptotic behavior of a Strang splitting scheme for damped second-order gradient systems' 
by Cristian Daniel Alecsa, Titus Pinta & Adrian Viorel

- the optimization algorithms can be found in the original paper
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
- a 2D sin type function that satisfies the PL inequality. It is similar to the function from the article 'A geometric integration approach to smooth optimisation: Foundations of the discrete gradient method' by M.J. Ehrhardt et. al. (we have put the coefficient 1/1.5 in front of the norm, in order to have plateaus around the unique global minimum. If one puts a value greater than 2 at the denominator then we would also have local minima around the global minimum)



The '.mw' files are made in Maple. These files contain the integrals from the RHS of the Strang splitting and Strang with Predictor-Corrector (for the 1D objective function) and the exact solution for the linear system of ODE's (in the 2D folder). Also, for Rosenbrock and 2D sin-type function, the '.mw' files contains the computed integrals that are used in '.py' files for the Newton methods for solving implicit systems.
