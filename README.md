# Passivity-FOs
Code associated with the paper "A Passivity Interpretation of Energy-Based ForcedOscillation Source Location Methods". This repository contains a folder ("Optimization_Functions") with the functions primarily used to solve the inference problems and perform some plotting. It also contains 5 files:

1) "Gen_Params.mat" - This is the set of generator parameters (machine + AVR) used in the testing.

2) "Load_Params.mat" - This is the set of load parameters (exponential + frequency dependent) used in the testing.

3) "NE39.m" - This file contains the PSAT static power flow data for the IEEE 39 bus power system

4) "Passivity_39.m" - This file is the main simulation file. It calls the system data, runs power flow via PSAT, computes the actual system response in the frequnecy domain, implemenets the inference testing, plots the dissipating power flow in the system, and computes current predictions.

5) "YaF_3rd_AVR.m" - This file builds the frequnecy response function associated with the generators.
