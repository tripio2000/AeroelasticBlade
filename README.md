# AeroelasticBlade
Aeroelastic simulation of a large horizontal wind turbine blade programmed in GNU Octave

## Description

Aeroelastic simulation of a single wind turbine blade under uniform flow and operating at a fixed angular velocity.

The aerodynamics is modelled with the Unsteady, three dimensional, vortex lattice method (UVLM3D).

The structural dyanmics is modelled using Rayleigh-Ritz approximation.

The developed model is based on a co-simulation of two subsystems that represent the aerodynamics and the structural dynamics of the blade.

The framework of the proposed simulations is highly modular so that the individual components can be replaced without modifying the overall organization. The proposed methodology is not restricted to periodic movements or linear movement equations. 

The simulations are developed in the time domain, therefore they allow the modeling of aeroelastic instabilities such as “flutter” and divergence.

This work is part of my Phd Thesis entitled "Estudio de palas inteligentes para mejorar la performance y la vida ́util de turbinas éolicas de gran
potencia y de eje horizontal"

The complete model is programmed in GNU Octave. Most of the operations are optimised for SIMD and the velocity induction is calculated using the "ocl" package for GPU acceleration.

## Contents

The main file is Uaero3D_GPU.m

Several parameters must be defined in the main file (TO DO change to external file for clarity)

### Input

The model requires descriptions for the aerodynamics and the structural dynamics of the blade.

In the folder INPUT the necessary files are included for testing.

### Output


## TO DO

Translate entirely to Englich language.



