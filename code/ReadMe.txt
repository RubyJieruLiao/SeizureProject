1.Main Pipeline: Simulation, Coupling Matrix Inference, Synchronizability, and Stability Results

(1) mainSingle.m ---- Runs the main pipeline for a signal parameter scenario.
(2) mainLoop.m ---- Runs the main pipeline for multiple parameter scenarios.

2.Seizure Simulation and Coupling Matrix

(1) oneDepileptor.m ---- Contains the equations for the 1D-Epileptor model.
(2) CouplingMatrix.m ---- Infers the weighted update matrix from the 1D-Epileptor model.

3.Calculate Synchronizability

(1) CovarianceUGaussianNet.m ---- Checks the necessary conditions for C.
(2) contCon2CovProjected.m ---- Calculates UcovarianceU.
(3) synchronizability.m ---- Computes /sigma^2 based on UcovarianceU.

4.Calcaulate Stability

(1) con2cov.m --- Calculatew covariance from /ncomp_tools
(2) Stability.m --- Computes /sigma^2_st based on covariance

Additional Files for Synchronizability and Stability
(1) maxrelerr.m
(2) negligible.m

5.Additional Code Files for Testing

(1) StabilityApproximation.m --- Calculates low approximations for stability results; not a function and can't be used as a .mat file
(2) PlotNodeFun.m --- Plots each node's values, such as z, BeAwayStability, and PushAwayStability.
(3) classify_nodes.m --- Identifies source and link nodes; algorithm based on the paper "Source-sink connectivity: a novel interictal EEG marker for seizure localization."