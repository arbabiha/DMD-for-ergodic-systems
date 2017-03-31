# DMD-for-ergodic-systems
A new variation of Dynamic Mode Decomposition for computing the Koopman eigenvalues for ergodic systems.
This folder contains the data and computation files used in
"Ergodic Theory, Dynamic Mode Decomposition and Computation of Koopman scpectral properties"
by Hassan Arbabi & Igor Mezic, 2016


The examples in the paper:

Lorenz_POD: computation of a POD basis for observables on chaotic Lorenz attractor,
PeriodicCavityFlow: computation of Koopman eigenvalues for periodic nonlinear flows using Hankel-DMD,
QPeriodicCavityFlow: computation of Koopman eigenvalues for periodic nonlinear flows using Exact Hankel-DMD,
VanDerPol_phase: computation of asymptotic phase for trajectories of Van der Pol oscillator.




The DMD routines in the +DMD folder:

algorithm suggested in the above paper= Hankel-DMD,
companion-matrix DMD = DMD_Vanilla,
SVD-enhanced DMD = DMD_Schmid,
Exact DMD	= EXACTDMD.

