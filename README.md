# P-ONE NuTau Project

Code for developing an algorithm to identify tau neutrinos using the P-ONE geometry. The scripts heavily rely on the IceCube software. To run the code in this repository the IceCube software should be installed. Most of the code was run on condor in Illume.

Steps to install the IceCube Software:

## Analysis

Contains code for analysis of the simulations.

- DOMandAngularAcceptance.py - Generate angular and wavelength acceptance of DOM simulated.
- ExpectedEventsInAYear.ipynb - Calculating the number of NuTau events expected to be detected in a year in P-ONE geometry, applying weights generated in NuGen
- I3PhotonDistributions.py - Plotting Wavelength, Positions, Azimuth and Zenith distributions(with and without the photon weight) of I3Photons simulated in CLSim.
- I3PhotonTimeDistibution.py - Plotting the time difference of the first and last photon to reach the DOM region in each event.
- NoiseTests.py - Plot the noise and physics hits distribution in each DOM in a single event.
- NuGenDistributions.py - Plotting simulated neutrino distributions of azimuth, zenith and energy distributions.
- PROPOSALAnalysis.py - Testing the working of PROPOSAL for Tau lepton.
- analysisHelpers.py - Helper functions to plot the geometry of the simulated detector and plot histograms(1D and 2D) of I3Photon scatter and direction distributions.

## BiGauss
This folder contains scripts to the initial analysis done in fitting a BiGauss and double BiGauss to each DOM in an event.   

- vertexDistances.py - Plotting distribution of distance between DOM and the tau creation and decay vertices.
- vertexDistances.sh, vertexDistances.submit - Scripts for running vertexDistances.py on Illume

## DOM

Angular acceptance and Wavelength acceptance data of the IceCube DOM and mDOM to be fed into the simulation.

## Flasher

Codes that simulate an isotropic flasher in icetray.

- genIsotropicFlashesSigma.py - Defining the flasher wavelengths and directions of flashes here.
- runCustomFlashes.py - Simulating multiple flashers to generate isotropic flashes.
- makePhotonsFromFlashes.py - Flashes from the flasher are converted to I3Photons.
- analysis.sh, analysis.submit - Scripts for running runCustomFlashes.py, makePhotonsFromFlashes.py  on Illume

## GCD

- PONE_Phase1.py produces .i3.gz file that replicates the phase 1 geometry of the P-ONE project.
- MediumCheck_GCD.py produces .i3.gz file that is used to test the flasher in the simulation.

## Graphing

- plots.py - Scatter, corner and histogram plots to plot the distributions of NuTau and background. The parameters used Time Difference, Width Difference, Amplitude Ratio and skew/k difference and log likelihood difference.

- plots2.py - Scatter plots for correlation of 2DeltaLLH with Amplitude asymmetry, energy of the lepton, brightness of events and tau length. The scatterPerEvent function considers the DOM with the largest 2DeltaLLH from a single event.

## IMINUIT

The analysis uses IMINUIT minimizer to find the best fit parameters of the biGauss/expGauss function.

- Debug.ipynb - Given the file numbers, frame numbers, string and DOM number the output from the minimizer is printed for that particular DOM. Histogram of the hits, along with the both a single expGuass fit and a double expGauss fit is plotted to verify how effective the parameters selected by the minimizer are.

- curveFit.py - Script that fits single and double biGauss/expGauss function to the hits distribution of each DOM and returns the best fit parameters. The list of parameters returned.
  - Double expGauss/BiGauss
    - log likelihood value
    - position 1
    - width 1
    - skew/k 1
    - amplitude 1
    - position 2
    - width 2
    - skew/k 2
    - amplitude 2

  - Single expGauss/biGauss
    - log likelihood value
    - position
    - width
    - skew/k
    - amplitude

- likelihoodHelpers.py - Functions for single and double exponentially modified gaussian/bifurcated gaussian. Contains also the functions for log likelihood functions.

- run_curveFit.py - Runs curveFit.py

- run_curveFit.submit, run_curveFit1.sh - Illume submission scripts to generate files using run_CurveFit.py

- Tau_E_analysis(update 1-4).ipynb - ipython notebooks that analyse the how well the separation between signal and background is. Single and double bifurcated gaussian are fit to each DOM in the event. Updates were made to the curveFit.py

- expGauss.ipynb - Identifying signal from background using single and double expGauss fits.

- expGauss_cuts.ipynb - Updated the DOM selection method in curveFit.py to include more number of NuTau events in the analysis. Separating NuTau events from background using the improved DOM selection.

- Weights.py - Functions for calculating the weight of the event.

- exitStatus.ipynb - Checking which last point in the curveFit.py the DOM has reached. The meaning of the assigned values is the following:
  - 0 - DOM has hits less than the assigned hit threshold
  - 1 - Less than 10 hits in the selected time window while calculating the mean
  - 2 - Bins in the histogram is less than a assigned bin threshold.
  - 3 - DOM made through the cuts previously imposed.

## Medium

- makeMediumSTRAW.py - generates the parameters to simulate water. The scatter and absorption values at different wavelengths are taken from Andreas Gaertner's analysis being done using STRAW data. The <cos(theta)> is calculated using values taken from Antares paper and eta from Matthew Man's analysis of STRAW data. The files generated are stored in STRAW_Andy_20200328_MattewEta

- STRAW_Andy_20200328_MattewEta is the folder containing properties of the medium to be fed into the simulation.

## Noise

Testing the noise injection in the simulation.

## NuTau_Simulation

Simulation chain for the P-ONE project in IceCube framework.

Produce Neutrinos ------> Propagate Secondary Particles ------> Geneterate Photons ------> Make Hits ------> Add Noise ------> Generate RecoPulses

The first three steps use code directly from the IceCube software. The last three steps are specifically written for P-ONE.

- step_1_neutrino_generator.py - Simulate neutrinos and propagate them using neutrino generator. Propagate the secondary particles produced in neutrino interactions in the medium using PROPOSAL.
- step_2_clsim.py - Uses CLSim, which produces and propagates photons; stores information about the photons that reach the DOM area.

Code written for P-ONE
- step_3_makeHits_pone.py - This takes information about the DOM( DOM acceptance and angular acceptance) used and area of the DOM to calculate the probability of registering a hit given a photon. Used the properties of IceCube mDOM. The properties of the mDOM are defined in MCPEConverter.py
- step_4_addNoise.py - Inject noise using the I3Module NoiseGenerator.py . Took the noise from the STRAW data and injected it directly into the simulation.
- step_5_makeRecoPulses.py - Here in the final step photoelectrons are converted into recopulses. Here hits within 3ns are merged. RecoPulseGenerator.py is used to generate the pulses.



## definingDOMs

Scripts used to find the double peak signature using the center of gravity.

## TimeResiduals

Codes for plotting the time residuals to verify the properties of the water in simulation against the properties of cascadia basin.
