import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_prominences


def generate_coincidence_data():
    ##------------------------------------------------------------------------
    ## Step 1. read in the data and separate it between channels
    ## Note: this is if CoMPASS was set to one file for all channels
    ## This script will need to be edited a little if each channel had its
    ## own .csv file
    ##------------------------------------------------------------------------
    path = "2022_04_01_AFITB10JJM_AmBe_coincidence"

    df = pd.read_csv(path+"/UNFILTERED/AmBe.csv" , sep = ';')

    ch0 = [] ; ch2 = [] ; ch4 = []
    for i in range(len(df["CHANNEL"])):
        if df["CHANNEL"][i] == 0:
            ch0.append([df["TIMETAG"][i],df["ENERGY"][i],df["ENERGYSHORT"][i]])
        elif df["CHANNEL"][i] == 2:
            ch2.append([df["TIMETAG"][i],df["ENERGY"][i],df["ENERGYSHORT"][i]])
        else:
            ch4.append([df["TIMETAG"][i],df["ENERGY"][i],df["ENERGYSHORT"][i]])


    np.save("ch0_data",ch0) ; np.save("ch2_data",ch2) ; np.save("ch4_data",ch4)

    ch0 = np.load("ch0_data.npy") ; ch2 = np.load("ch2_data.npy") ; ch4 = np.load("ch4_data.npy")
    ##------------------------------------------------------------------------
    ## Step 2. Isolate the 480 keV gamma peak that is present in the LaBr3
    ##         spectrum for each channel
    ##------------------------------------------------------------------------
    data,bins = np.histogram(ch2[:,1],500)
    peaks, _ = find_peaks(data) # locates the peaks in the binned spectrum
    prominences = peak_prominences(data, peaks)[0]
    maxprominences1 = np.argmax(prominences) # pulls out the index of the 1st peak
    prominences[maxprominences1] = 0
    maxprominences2 = np.argmax(prominences) # pulls out the index of the 2nd peak
    ## located the centroid of the 480 keV full energy peak
    centroid2 = peaks[maxprominences2] # 1st peak is low channel noise so we want 2nd
    centroid2 = bins[centroid2]
    ##------------------------------------------------------------------------
    ## Ch4 begins
    ##------------------------------------------------------------------------
    data,bins = np.histogram(ch4[:,1],500)
    peaks, _ = find_peaks(data) # locates the peaks in the binned spectrum
    prominences = peak_prominences(data, peaks)[0]
    maxprominences1 = np.argmax(prominences) # pulls out the index of the 1st peak
    prominences[maxprominences1] = 0
    maxprominences2 = np.argmax(prominences) # pulls out the index of the 2nd peak
    # Ch4 full energy peak acquired
    centroid4 = peaks[maxprominences2] # 1st peak is low channel noise so we want 2nd
    centroid4 = bins[centroid4]
    ##------------------------------------------------------------------------
    ## Step 3. apply a cut to Ch2/4 data sets to only include 480 keV events
    ##------------------------------------------------------------------------
    cut2, stop2 = centroid2-120, centroid2+120 # +/- value depends on detector resolution
    cut4, stop4 = centroid4-120, centroid4+120

    ch2 = np.array([x for x in ch2 if cut2 <= x[1] <= stop2])
    ch4 = np.array([x for x in ch4 if cut4 <= x[1] <= stop4])
    ##------------------------------------------------------------------------
    ## Step 4. Now check the time of each channel's events and match up
    ##         coincidence between boron-loaded (CH0) and either LaBr3
    ##------------------------------------------------------------------------
    min_time = min(min(ch2[:,0]),min(ch4[:,0]))
    max_time = max(max(ch2[:,0]),max(ch4[:,0]))

    ch0 = np.array([x for x in ch0 if min_time <= x[0] <= max_time])

    # convert times to minutes

    time0 = (ch0[:,0]/1e12)/60
    time2 = (ch2[:,0]/1e12)/60
    time4 = (ch4[:,0]/1e12)/60

    time0 = np.array([round(x,5) for x in time0])
    time2 = np.array([round(x,5) for x in time2])
    time4 = np.array([round(x,5) for x in time4])

    ch02 = np.array([ch0[i] for i in range(len(time0)) if time0[i] in time2])
    ch04 = np.array([ch0[i] for i in range(len(time0)) if time0[i] in time4])



    np.save("ch02",ch02) ; np.save("ch04",ch04)
