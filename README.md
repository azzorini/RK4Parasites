# RK4Parasites

## Introduction

Simulations of parasite produced marine diseases by using fourth order Runge-Kutta method. This simulations are all based on the paper [Modelling parasite-produced marine diseases: The case of the mass mortality event of Pinna nobilis](https://www.sciencedirect.com/science/article/pii/S030438002100260X). In this paper basically a four dimmensional dynamcial system is proposed and then some reductions and approximations are applied.

## Files

In this repository yoou can find the following files:

1. *SIRP_exact.cpp:* In this file the complete four order dynamical system is simulated.
2. *SIRP_exact_red.cpp:* As it is shown in the reference paper this system has a conserved quantity. This quantity is used to reduce the dimensionality of the system in an exact way. This code simulates this exactly reduced system.
3. *SIRP_approx_red.cpp:* In the paper an approximation is introduced to make the exact reduced system more convenient. This code simulate this approximation to the exact reduction.
4. *SIRP_fastslow_approx:* Eventually a time scale separation was done. This code simulate this approximation.
