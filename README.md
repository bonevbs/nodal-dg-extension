# nodal-dg-extenstion
A simple extension for the nodal-dg code by Jan S. Hesthaven and Tim Warburton

## Description
This library includes extensions to the original nodal-dg code by Jan S. Hesthaven and Tim Warburton. So far, these extensions include:
* High-order continuous Galerkin (CG/FEM) method on triangular meshes using the nodal-dg datastructures
* DG method for the high-contrast Poisson problem
* IPDG method for linear elasticity in 2D
* Test routines designed to check behaviour under h- and p-refinement

## Setup
The code needs the original nodal-dg code which can be found here: [nodal-dg](https://github.com/tcew/nodal-dg)
To run the codes, simply include both the original nodal-dg folder and the nodal-dg-extension folder in the Matlab path.
