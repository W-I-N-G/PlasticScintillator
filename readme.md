# Plastic Scintillator Work

**to use the automated waveform analysis script, "initialize_routine.py", please
ensure you have the addcopyfighandler module installed,
https://pypi.org/project/addcopyfighandler/

We ask that you please cite our paper, NIM link, if you used the automated script

Thank you!

The automated routine supports Na22, Cs137, and AmBe dataset analysis. There are
several outputs from the process,
----------------------------------------------------
Plot Outputs (saved in the plots/ dir)
----------------------------------------------------
1. Na22/Cs137/AmBe PSD plot
2. Na22/Cs137 uncalibrated spectrum
3. Na22/Cs137/AmBe waveform plot with gates labeled
4. AmBe Figure of Merit plot
----------------------------------------------------
Text file Outputs (saved in your data/outputs/ dir)
----------------------------------------------------
1. dataset_source.txt
      - this file will hold the light calibration data in units keVee
----------------------------------------------------

First, ensure the following directory exists

1. data/

in the data/ directory should be the CoMPASS data you'd like to analyze. Please
make sure that the different datasets have the source located in the file name.
For example,

- data/ambe_data.filetype
- data/02March2022_ambe_waveforms.filetype
- data/Na-22-gammadata.filetype

This is important because the algorithm will perform a keyword search to capture
the correct dataset for each source. To speed things up, the algorithm will save
the CoMPASS data as .NPY file types.

To execute the process simply call the "initialize_routine.py" script.
This script will automatically generate the spectra for each source,
locate the Compton Edge(s), fit a Gaussian to the(those) edge(s), and calibrate
accordingly. For AmBe data, it will also generate the pulse-shape discrimination
(PSD) plot. This plot is useful for separating Gamma and Neutron events.
A PSD cut is applied on the for Gamma/Neutron events on the interval
[450, 550] keV and a double-Gaussian function is applied as a fit and the
Figure of Merit (FoM) is maximized as a function of the PSD short gate value.
As the process is working several values of interest will be printed to the
user's terminal. An example of executing the process has been provided below.
