# navier-stokes-rust

## Introduction

This project is a 2D fluid dynamics simulation built from the ground up in Rust. It solves the incompressible Navier-Stokes equations on a grid to simulate the motion of viscous fluids. It is part of my attempt to learn the Rust programming language.

## How to use
Make sure you have Rust installed properly on your local machine. This can be done by visiting the installation instructions here:

https://doc.rust-lang.org/book/ch01-01-installation.html

Once Rust is installed locally, clone the repo locally and build/run it with:

`cargo run --release` 

## Mathematical Foundation
The simulation uses a classic computational fluid dynamics (CFD) approach. I used the incompressible Navier-Stokes equations which govern the relationship between a fluid's velocity, pressure, density, and viscosity. Since they are a set of complex partial differential equations, we can't really solve for them directly -- we can just use computational methods to find approximate solutions. 

Since I focused on imcompressible fluids (density is constant, like water), I only needed to worry about two main equations:

#### Momentum Equation
This is basically Newton's Second Law ($F = ma$) rewritten for a fluid. It describes how the velocity of the fluid changes over time due to forces:

$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2 \mathbf{u} + \mathbf{g}$

#### Continuity Equation (Incompressibilty)
This just enforces the conservation of mass. It says that the *divergence* of the velocity field is zero. For any infinitesimally small volume within the fluid, the amount of fluid flowing in must equal the amount flowing out. The fluid can't be compressed (incompressible), created, or destroyed in any way.

$\nabla \cdot \mathbf{u} = 0$

## Bridging to the computer
We have to do a bit of "bridging" from the continuous ways of calculus to the discrete ways of computers. Basically, I just converted the derivatives in the equations into a large system of linear equations (finite difference formulas) using something called the Finite Difference Method. I'll post some sources to read about it at the end of this document. 

## The Grid
For the grid of cells, I decided to use a staggered (MAC) grid instead of a collocated grid. Pressure, $p$, is stored at the center of the cell, while the horizontal and vertical velocities, $u$ and $v$ respectively, are stored on the vertical and horizontal faces, respectively. 

This does add complexity, sure, but it gives us an advantage! When we calculate the finite difference, the variables are exactly where we need them:
- Pressure Force --> To find the pressure force acting on the horiztonal velocity, $u$, we need the pressure difference between the cell to its left and the cell to its right. On our staggered grid, the $u$ component lives exactly between the two pressure points. 
- Divergence --> To enforce the incompressibility condition ($\nabla \cdot \mathbf{u} = 0$) at the center of the cell, we need to know the velocity flowing in and out through its faces. Our staggered grid places the velocity components right on these faces.

## Solver Algorithm
We use the projection method as the engine of the simulation. The diction for this can be a bit daunting those not familiar with the mathematics, but I'll put some sources at the end to read up on it.

It is a guess-and-correct strategy for our algorithm. Here it is broken down into three steps:
1. Advection --> A Semi-Lagrangian scheme is used to advect the velocity field, to movie it along with the flow. Bilinear interpolation was used to find the values between grid points. 
2. Pressure Solve --> The pressure field is calculated by solving a Poisson equation using the Jacobi iterative method. We do this to ensure the velocity field remains divergence-free and correct throughout the simulation. This ties in with conserving the mass throughout the entirety of our simulation.

Poisson Equation:

$$\nabla^2 p^{n+1} = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*$$

3. Projection --> The final velocity field is corrected by subtracting the pressure gradient. 

Correction equation:

$$\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla p^{n+1}$$

That forms the entirety of our loop. We start with a valid velocity, make a guess based on advection, solve for the pressure that keeps it divergence-free and enforces mass conservation, and then use that pressure to correct our initial valid velocity. We do this for each time step. 


## Implementing in Rust

My implementation in Rust centered around a `FluidGrid` struct to hold all the simulation data. I organized the actual steps by just following the mathematical flow dictated above!

To display our simulation, I used the `piston_window` library for real-time rendering. 

Currently, the only visualization mode is a vector field plot that renders velocity vectors as arrows to show the direction and magnitude of the flow. I would like to implement a color gradient to display the speed of the flow as well.

Was Rust the right language for this? Maybe not. Was it a fun way to dive into the intricacies of Rust? Yes!


## Sources

This textbook I picked up was quite helpful:
https://www.amazon.com/Simulation-Computer-Graphics-Robert-Bridson/dp/1482232839

and this paper was a good read:

https://cg.informatik.uni-freiburg.de/intern/seminar/gridFluids_fluid-EulerParticle.pdf



