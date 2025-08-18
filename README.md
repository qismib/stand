# stand
Simulation and Trotter Analysis of Neutrino Dynamics

## /commutators
This folder contains scripts that compute (nested) commutators relevant to my analysis of Trotter error.\\
The code makes large use of `CommutatorTool.py` class.

## /simulations
This folder contains scripts that simulate the dynamics of quantum systems using `qutip`.

## TODO
- implement CommutatorTool 2.0 based on `numpy`.
- numerical estimate of first and second order additive errors, comparing exact propagator against the Trotterized version and comparison to analytical bound.
- ...