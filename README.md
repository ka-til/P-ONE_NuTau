# P-ONE NuTau Project

Note - Branch BiGauss has the latest updates, yet tp be merged with the master

## Analysis

Contains code for analysis of the simulations.

## DOM

Contains Angular acceptance and Wavelength acceptance of the IceCube DOM to be fed into the simulation.

## Flasher

Codes that simulate an isotropic flasher in icetray.

## GCD

PONE_Phase1.py produces .i3.gz file that replicates the phase 1 geometry of the P-ONE project.
MediumCheck_GCD.py produces .i3.gz file that is used to test the flasher in the simulation.

## Medium

makeMediumSTRAW.py generates the parameters to simulate water.
STRAW_Andy_20200328_MattewEta is the folder containing properties of the medium to be fed into the simulation.

## Noise

Simulation and analysis files for noise. Move the files, once the code is finalized.

## NuTau_Simulation

Simulation chain for the P-ONE project in IceCube framework.

Produce Neutrinos ------> Propagate Secondary Particles ------> Geneterate Photons ------> Make Hits ------> Add Noise

## TimeResiduals

Codes for plotting the time residuals to verify the properties of the water in simulation against the properties of cascadia basin.

