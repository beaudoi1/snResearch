snResearch — Circumstellar Material in the UV Spectra of Type Ia Supernovae

Undergraduate research code from the University of Illinois Department of Astronomy (2014–2015), written by Aaron Beaudoin under the supervision of Prof. Ryan J. Foley.

The project measured narrow ultraviolet absorption lines in Hubble Space Telescope spectra of 8 Type Ia supernovae (39 spectra). The goal was to test whether gas around the exploding star changes over time. Earlier studies relied on the sodium (Na I) doublet. This work was a first look at Mg II, Mg I, and Fe II in the UV, which probe gas that Na I can't.

What the code does

spAnalyze.py (~2,000 lines, Python) is an interactive analysis tool built to run in IPython. For each spectrum it:

Fits the continuum with a 3rd- or 5th-order B-spline. The order is chosen by signal-to-noise, and the breakpoints can be added, removed, and saved by hand. Dividing by the fit gives a normalized spectrum.
Fits each absorption line with a Gaussian using least-squares minimization, then integrates the Gaussian to get the equivalent width (EW).
Gaussian fits are more robust than direct integration on low-S/N data, where the chosen wavelength window biases the result.
Separates overlapping lines. When the Milky Way and host-galaxy Mg II doublets overlap in wavelength, it fits both doublets at once. Fe II blends are handled the same way.
Estimates systematic error with a Monte Carlo simulation (500 iterations). Each iteration shifts the spline breakpoints, re-derives the continuum, and re-measures every line. The scatter becomes the systematic uncertainty, which is combined with the statistical error from the fit.
Sets 3σ detection limits for lines that aren't detected, based on local noise.
Checks fit quality with reduced chi-square tests.
Saves the results and produces output: parameters, breakpoints, and measurements go to per-epoch files; the tool also plots EW vs. epoch and exports LaTeX tables.

Lines measured (rest wavelengths, Å):

Ion	Lines
Fe II	2344.21, 2374.46, 2382.76, 2586.65, 2600.17
Mg II	2796.35, 2803.53
Mg I	2852.96
Sample

SN 1992A, SN 2011by, SN 2011ek, SN 2011fe, SN 2011iv, SN 2012cg, SN 2013dy, SN 2014J.

All were observed by HST: STIS for every object except SN 1992A, which used FOS. The sample includes every SN Ia with at least two high-resolution (R > 500), high-S/N spectra covering the Mg II doublet.

Preliminary findings (unfinished draft)
Mg II absorption was measurable in 5 of the 8 SNe.
The other three (SNe 1992A, 2011iv, 2011ek) showed only Milky Way lines. That is consistent with their early-type hosts or large offset from the host galaxy.
Every SN with detectable Mg II also showed Mg I λ2852.
Fe II absorption was detected only in SN 2011fe.

The paper was not completed, so treat these results as preliminary.

Repository layout
Path	Contents
spAnalyze.py	Main analysis tool: continuum fitting, Gaussian EW fits, Monte Carlo errors, limits, plotting
orgData.py, makeTables.py	Collect measurements across SNe and build LaTeX tables
limits.py, conv06X.py	Detection-limit tests; convolving spectra to match another instrument's resolution
plot.py, compPlot.py, resIon.py, mgivsii.py, contPlots.py	Figures: EW over time, change in EW, Mg I vs. Mg II, continuum fits
sn2011fe/, sn2013dy/	Early prototype scripts (B-spline fitting, Gaussian tests, Monte Carlo tests)
sn*/Data/	Per-supernova input spectra (.flm), breakpoints, fit parameters, EW results, and plots
paper.tex, data_tables.tex, obslog.tex	Draft paper and tables
Tech

Python 2.7 · NumPy · SciPy (interpolation, optimization, integration, statistics) · Astropy (tables, convolution) · Matplotlib · pyspeckit (mpfit) · LaTeX

Note: written in 2014–2015 for Python 2 (it uses print statements). Running it today requires Python 2.7 or porting to Python 3.

Usage

The tool is interactive, not a batch pipeline. From the repo root, start IPython and run:

python
import spAnalyze as a
a.pickSN()   # choose a supernova and load its parameters

From there, the functions listed in the header of spAnalyze.py step through the analysis: breakpoints, EW fits, Monte Carlo, limits, and plots.snResearch
==========
