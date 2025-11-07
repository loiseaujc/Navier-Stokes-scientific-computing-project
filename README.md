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
