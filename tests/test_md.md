# Contents
- [Contents](#contents)
  - [1 Error Estimate Formulas](#1-error-estimate-formulas)
  - [2 Symbolical Commutator Evaluators](#2-symbolical-commutator-evaluators)
    - [2.1 Heisenberg Hamiltonian Terms](#21-heisenberg-hamiltonian-terms)
      - [2.1.1 The `gen_Heisenberg_terms` function](#211-the-gen_heisenberg_terms-function)
    - [2.2 Rules for Commutators of Pauli Strings](#22-rules-for-commutators-of-pauli-strings)
    - [2.3 Evaluation of Commutators: the `comm_pstr` function](#23-evaluation-of-commutators-the-comm_pstr-function)
      - [2.3.1 $N=3$ tests](#231-n3-tests)
  - [3 First Order Single-Step Error](#3-first-order-single-step-error)
    - [3.1 Compute Error Terms Numerically](#31-compute-error-terms-numerically)
    - [3.2 $e\_1$ decomposition](#32-e_1-decomposition)
    - [3.3 Numerical Analysis of Error Scalings](#33-numerical-analysis-of-error-scalings)
      - [3.3.1 $dt$ scaling](#331-dt-scaling)
      - [3.3.2 $N$ scaling: there are caveats](#332-n-scaling-there-are-caveats)
    - [3.4 Peaked Neutrino Beam](#34-peaked-neutrino-beam)
      - [3.4.1 Visualization of neutrino beam](#341-visualization-of-neutrino-beam)
      - [3.4.2 $dt$ scaling](#342-dt-scaling)
      - [3.4.3 $N$ scaling: first approach (log-log)](#343-n-scaling-first-approach-log-log)
      - [3.4.4 $N$ scaling: second approach (alt)](#344-n-scaling-second-approach-alt)
      - [3.4.5 $N$ scaling: third approach (count terms to estimate $\\gamma'$) (TODO)](#345-n-scaling-third-approach-count-terms-to-estimate-gamma-todo)
      - [3.4.6 Compare $N$ scalings results](#346-compare-n-scalings-results)
  - [4 First Order Finite Time Scalings](#4-first-order-finite-time-scalings)
  - [5 Second Order Single-Step Error](#5-second-order-single-step-error)
    - [5.1 Compute Error Terms Numerically](#51-compute-error-terms-numerically)
    - [5.2 $e\_2$ decomposition](#52-e_2-decomposition)
    - [5.3 Numerical analysis of error scalings](#53-numerical-analysis-of-error-scalings)
      - [5.3.1 $dt$ scaling](#531-dt-scaling)
    - [5.4 Peaked neutrino beam](#54-peaked-neutrino-beam)
      - [5.4.1 $dt$ scaling](#541-dt-scaling)
      - [5.4.2 $N$ scaling: first approach (log-log)](#542-n-scaling-first-approach-log-log)
      - [5.4.3 $N$ scaling: second approach (alt)](#543-n-scaling-second-approach-alt)
      - [5.4.4 $N$ scaling: third approach (count terms to estimate $\\gamma'$) (TODO)](#544-n-scaling-third-approach-count-terms-to-estimate-gamma-todo)
      - [5.4.5 Compare $N$ scalings results](#545-compare-n-scalings-results)






## 1 Error Estimate Formulas



## 2 Symbolical Commutator Evaluators

### 2.1 Heisenberg Hamiltonian Terms

#### 2.1.1 The `gen_Heisenberg_terms` function

### 2.2 Rules for Commutators of Pauli Strings

### 2.3 Evaluation of Commutators: the `comm_pstr` function

#### 2.3.1 $N=3$ tests


## 3 First Order Single-Step Error

### 3.1 Compute Error Terms Numerically

### 3.2 $e_1$ decomposition

### 3.3 Numerical Analysis of Error Scalings

#### 3.3.1 $dt$ scaling

#### 3.3.2 $N$ scaling: there are caveats

### 3.4 Peaked Neutrino Beam

#### 3.4.1 Visualization of neutrino beam

#### 3.4.2 $dt$ scaling

#### 3.4.3 $N$ scaling: first approach (log-log)

#### 3.4.4 $N$ scaling: second approach (alt)

#### 3.4.5 $N$ scaling: third approach (count terms to estimate $\gamma'$) (TODO)

#### 3.4.6 Compare $N$ scalings results


## 4 First Order Finite Time Scalings


## 5 Second Order Single-Step Error

### 5.1 Compute Error Terms Numerically

### 5.2 $e_2$ decomposition

### 5.3 Numerical analysis of error scalings

#### 5.3.1 $dt$ scaling

### 5.4 Peaked neutrino beam

#### 5.4.1 $dt$ scaling

#### 5.4.2 $N$ scaling: first approach (log-log)

#### 5.4.3 $N$ scaling: second approach (alt)

#### 5.4.4 $N$ scaling: third approach (count terms to estimate $\gamma'$) (TODO)

#### 5.4.5 Compare $N$ scalings results
