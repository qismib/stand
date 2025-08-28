# stand
Simulation and Trotter Analysis of Neutrino Dynamics

## /commutators
This folder contains scripts that compute (nested) commutators relevant to my analysis of Trotter error.\\
The code makes large use of `CommutatorTool.py` class.

## /simulations
This folder contains scripts that simulate the dynamics of quantum systems using `qutip`.

## TODO
- implement CommutatorTool 2.0 based on `numpy`.
- work out N=4 results in order to have means to compare results directly with https://arxiv.org/pdf/2207.03189.
- implement definition of coupling constants in terms of physical constants being careful to not miss any conversion factor. 
- ...