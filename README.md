# Scientific Computing Project

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)

| **Program** | **Title** |
| :---: | :--: |
| Master Factory of the Future | Introduction to CFD |

## Description

This repository contains the skeleton code for the scientific computing project *Introduction to Computational Fluid Dynamics* for the International Master Program *Factory of the Future* at [Arts & Métiers Institute of Technology](https://artsetmetiers.fr/en/mecanique-energie-et-ingenierie). Inspired by the [12 steps to Navier-Stokes](https://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/), it aims at introducing students to the basics of computational fluid dynamics. It also aims at familiarizing students with good software development practices including `git` (and GitHub), unit testing and CI/CD.

## Practical considerations

The project is to be completed using the `Fortran` programming language. While one might argue that `your-favorite-programming-language` is a better choice for whatever reasons, `Fortran` has a number of built-in features justifying this choice. These include:

- `Fortran` has been designed from the ground up for computationally intensive applications in science and engineering.
- It is a statically and strongly typed langauge, allowing the compiler to catch many programming errors early on for you. It also allows the compiler to generate extremely efficient binary code.
- `Fortran` is relatively small language that is surprisingly easy to learn and use.
- It is a natively parallel programming language with intuitive array-like syntax. It also supports both object-oriented and functional programming paradigms.
- Many CFD solvers actually use `Fortran` under the hood and this project let students get somewhat familiar with it.

More information about `Fortran` can be found online at [fortran-lang.org](https://fortran-lang.org/) and the dedicated [discourse](https://fortran-lang.discourse.group/).
Relying on the modern `Fortran` ecosystem, the project will make heavy use of the following tools and libraries:

- `fpm` : the `Fortran` package manager ([documentation](https://fpm.fortran-lang.org/)), available for all operating systems (Windows, macOS, Linux).
- `stdlib` : the `Fortran` standard library ([documentation](https://stdlib.fortran-lang.org/)) providing high-performance yet easy-to-use numerical linear algebra functionalities.
- `test-drive` : a simple Fortran package for unit-testing directly supported by the Fortran package manager.
- `git` : A version control system to keep track of the incremental changes made to the code base.

To install any of these tools, please refer to the appropriate documentation and whatever package manager available in your operating system.

## Organisation

As stated before, this project is inspired by the [12 steps to Navier-Stokes](https://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/). As such, it uses a similar progressive learning curve.

### One-dimensional problems

As a starting point, students will familiarize themselves with both `Fortran` and standard numerical methods on simple one-dimensional problems.

**Step 1 - Steady heat equation.** The `Hello World!` of partial differential equations, it reads

$$
-\dfrac{d^2 u}{d x^2} = f(x).
$$

This equation will be discretized with the standard second-order accurate central finite-difference scheme on the unit segment $L = \left[ 0, 1 \right\]$ and uniform grid spacing leading to a tridiagonal matrix representation of this differential operator. The corresponding algebraic system of equations will be solved with a standard dense solver and the [Thomas algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm).

*Learning objectives.*

- Discretization of differential operators
    - Derivation of the finite-difference schemes.
    - Relation between the order of finite-difference scheme and the discretization error.
- Numerical linear agebra
    - Performances of dense vs. structure-aware linear solvers.

**Step 2 - Unsteady heat equation.** Having solved for the steady solution, the next natural step is to simulate the transients. For that purpose, we'll consider the unsteady heat equation

$$
\dfrac{\partial u}{\partial t} = \dfrac{\partial^2 u}{\partial x^2} + f.
$$

The spatial differential operator will be discretized as before. The focus is on the time-integration schemes. Three such schemes will be considered:

- First order explicit Euler
- First order implicit Euler
- Second order Crank-Nicholson

Both the implicit Euler and the Crank-Nicholson implementations will leverage the knowledge acquired during *step 1* on linear solvers.

*Learning objectives.*

- Time-integration of partial differential equations
    - Courant-Friedrich-Levy condition number
    - Numerical stability of an integration scheme
    - Pros and cons of explicit vs implicit temporal schemes.

**Step 3 - Advection-diffusion equation.**

**Step 4 - Nonlinear Burger's equation.**

### Two-dimensional problems

**Step 1 - Poisson's equation.**

