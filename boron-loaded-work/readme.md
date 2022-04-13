# Boron Loaded Work

Here are the post-processing scripts and data files used with the boron-loaded
plastic scintillator. The coincidence data set was too large to include on GitHub,
but it available upon request.

Directories:

1. 2022_04_01_AFITB10JJM_AmBe_BisBox/ -- AmBe data that was recorded while the
   plastic scintillator was enclosed in a bismuth box

2. 2022_04_04_AFITB10JJM_Am241/ -- Am241 data that was recorded for calibration

Python scripts:

1. gamma_data.py -- this script analyzes the Am241 spectrum, isolates the Am241
   59.5 keV photopeak, and applies a Gaussian + exponential background best-fit
   it will generate a plot of the entire Am241 spectrum as well as a zoomed in
   plot around the 59.5 keV photopeak with the best-fit plotted over the signal

2. ambe_data.py -- this script analyzes the AmBe_bisbox data and generates a
   PSD vs. Light histogram which showcases the thermal neutron response of the
   detector. It also conducts the coincidence data analysis that was generated
   using 2 LaBr3 detectors alongside the boron-loaded detector. This script calls
   on the coincidence.py script. The outputs of this script are light yield data,
   a PSD vs. Light histogram plot, a Time vs. Light histogram plot, a plot of the
   coincidence signal without timing considerations, and a plot of the neutron
   centroid feature location versus time.

3. coincidence.py -- here is the script responsible for scanning the coincidence
   data set and lining up any events that were coincident between the boron-loaded
   scintillator and either one of the LaBr3 scintillators. It will generate five
   .NPY data files,

   3.1. ch0.npy -- a Nx3 array ("TIME", "ENERGY", "ENERGYSHORT") for events in ch0
   3.2. ch2.npy -- same as 3.1 for events in ch2
   3.3. ch4.npy -- same as 3.1 for events in ch4
   3.4. ch02.npy -- coincidence events between ch0 and ch2 Nx3 array ("TIME", "ENERGY", "ENERGYSHORT")
   3.5. ch04.npy -- same as 3.4 for events between ch0 and ch4

   ambe_data.py will read in the ch02.npy and ch04.npy files for further analysis

4. gaussian.py -- holds the various gaussian fit functions that were used or sampled
   throughout this work
