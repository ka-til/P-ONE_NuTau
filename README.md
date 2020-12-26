# P-ONE NuTau Project

Code for developing an algorithm to identify tau neutrinos using the P-ONE geometry. The scripts heavily rely on the IceCube software. To run the code in this repository the IceCube software should be installed. Most of the code was run on condor in Illume.

Steps to install the IceCube Software:

## Analysis

Contains code for analysis of the simulations.

DOMandAngularAcceptance.py - Generate angular and wavelength acceptance of DOM simulated.
ExpectedEventsInAYear.ipynb - Calculating the number of NuTau events expected to be detected in a year in P-ONE geometry, applying weights generated in NuGen
I3PhotonDistributions.py - Plotting Wavelength, Positions, Azimuth and Zenith distributions(with and without the photon weight) of I3Photons simulated in CLSim.
I3PhotonTimeDistibution.py - Plotting the time difference of the first and last photon to reach the DOM region in each event.
NoiseTests.py - Plot the noise and physics hits distribution in each DOM in a single event.
NuGenDistributions.py - Plotting simulated neutrino distributions of azimuth, zenith and energy distributions.
PROPOSALAnalysis.py - Testing the working of PROPOSAL for Tau lepton.
analysisHelpers.py - Helper functions to plot the geometry of the simulated detector and plot histograms(1D and 2D) of I3Photon scatter and direction distributions.

## BiGauss
This folder contains scripts to the initial analysis done in fitting a BiGauss and double BiGauss to each DOM in an event.   

vertexDistances.py - Plotting distribution of distance between DOM and the tau creation and decay vertices.
vertexDistances.sh, vertexDistances.submit - Scripts for running vertexDistances.py on Illume

## DOM

Angular acceptance and Wavelength acceptance data of the IceCube DOM and mDOM to be fed into the simulation.

## Flasher

Codes that simulate an isotropic flasher in icetray.

genIsotropicFlashesSigma.py - Defining the flasher wavelengths and directions of flashes here.
runCustomFlashes.py - Simulating multiple flashers to generate isotropic flashes.
makePhotonsFromFlashes.py - Flashes from the flasher are converted to I3Photons.
analysis.sh, analysis.submit - Scripts for running runCustomFlashes.py, makePhotonsFromFlashes.py  on Illume

## GCD

PONE_Phase1.py produces .i3.gz file that replicates the phase 1 geometry of the P-ONE project.
MediumCheck_GCD.py produces .i3.gz file that is used to test the flasher in the simulation.

## Graphing

plots.py - Scatter, corner and histogram plots to plot the distributions of NuTau and background. The parameters used Time Difference, Width Difference, Amplitude Ratio and skew/k difference and log likelihood difference.

plots2.py - Scatter plots for correlation of 2DeltaLLH with Amplitude asymmetry, energy of the lepton, brightness of events and tau length. The scatterPerEvent function considers the DOM with the largest 2DeltaLLH from a single event.

## IMINUIT

The analysis uses IMINUIT minimizer to find the best fit parameters of the biGauss/expGauss function.

Debug.ipynb - Given the file numbers, frame numbers, string and DOM number the output from the minimizer is printed for that particular DOM. Histogram of the hits, along with the both a single expGuass fit and a double expGauss fit is plotted to verify how effective the parameters selected by the minimizer are.

curveFit.py - Script that fits single and double biGauss/expGauss function to the hits distribution of each DOM and returns the best fit parameters.

run_curveFit.py - Runs curveFit.py

run_curveFit.submit, run_curveFit1.sh - Illume submission scripts to generate files using run_CurveFit.py

Tau_E_analysis(update 1-4).ipynb - ipython notebooks that analyse the how well the separation between signal and background is. Single and double bifurcated gaussian are fit to each DOM in the event. Updates were made to the curveFit.py

expGauss.ipynb - Identifying signal from background using single and double expGauss fits.

expGauss_cuts.ipynb - Updated the DOM selection method in curveFit.py to include more number of NuTau events in the analysis. Separating NuTau events from background using the improved DOM selection.

Weight

## Medium

makeMediumSTRAW.py generates the parameters to simulate water.
STRAW_Andy_20200328_MattewEta is the folder containing properties of the medium to be fed into the simulation.

## Noise

Simulation and analysis files for noise. Move the files, once the code is finalized.

## NuTau_Simulation

Simulation chain for the P-ONE project in IceCube framework.

Produce Neutrinos ------> Propagate Secondary Particles ------> Geneterate Photons ------> Make Hits ------> Add Noise ------> Generate RecoPulses

## TimeResiduals

Codes for plotting the time residuals to verify the properties of the water in simulation against the properties of cascadia basin.
