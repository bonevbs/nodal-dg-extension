# nodal-dg-extension
A simple extension for the discontinuous Galerkin code [nodal-dg](https://github.com/tcew/nodal-dg) by Jan S. Hesthaven and Tim Warburton.

## Description
This library includes extensions to the original nodal-dg code. So far, these extensions include:
* High-order continuous Galerkin (CG/FEM) method on triangular meshes using the nodal-dg datastructures
* DG method for the high-contrast Poisson problem
* IPDG method for linear elasticity in 2D
* Routines for generating a nested dissection for the structured elimination of the resulting matrices
* Plotting routines designed to visualize the results
* Test routines designed to check behaviour under h- and p-refinement

## Setup
The code needs the original nodal-dg code which can be found here: [nodal-dg](https://github.com/tcew/nodal-dg)
To run the codes, simply include both the original nodal-dg folder and the nodal-dg-extension folder in the Matlab path.

## References
Hesthaven, J. S., & Warburton, T. (2008). Nodal Discontinuous Galerkin Methods. https://doi.org/10.1007/978-0-387-72067-8
