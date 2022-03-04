# Import packages
import os
import numpy as np
import matplotlib.pyplot as plt
import addcopyfighandler
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D
from scipy.signal import find_peaks, peak_prominences
import matplotlib as mpl


#%%
def PlotWaveforms(waveforms,numberofwaveforms,nSamplesBaseline,nSamplesShort,nSamplesLong):
    # Plot the first few waveforms
    plt.figure(constrained_layout=True)
    for i in range(numberofwaveforms):
        plt.plot(waveforms[i,:])

    # Plot the pregate, short gate and long gate
    samples = len(waveforms[0,:])
    baseline_plot = np.ones(samples)*15000
    long_plot = np.ones(samples)*15000
    short_plot = np.ones(samples)*15000
    baseline_plot[:nSamplesBaseline] = 15500
    short_plot[nSamplesBaseline:nSamplesBaseline+nSamplesShort] = 16000
    long_plot[nSamplesBaseline:nSamplesBaseline+nSamplesLong] = 16500
    plt.plot(baseline_plot,label='Baseline Gate',color='b')
    plt.plot(short_plot,label='Short Gate',color='g')
    plt.plot(long_plot,label='Long Gate',color='r')
    plt.xlabel('Sample')
    plt.ylabel('ADC')
    plt.legend()

#%%
def CalculateElongEshortPSD(waveforms,nSamplesBaseline,nSamplesShort,nSamplesLong):
    baseline = np.mean(waveforms[:,:nSamplesBaseline],axis=1)
    eshort = np.sum(baseline[:,None] - waveforms[:,nSamplesBaseline:nSamplesBaseline+nSamplesShort],axis=1)
    elong = np.sum(baseline[:,None] - waveforms[:,nSamplesBaseline:nSamplesBaseline+nSamplesLong],axis=1)
    PSD = (elong-eshort)/elong
    data = np.hstack((elong.reshape(len(elong),1),eshort.reshape(len(eshort),1),PSD.reshape(len(PSD),1)))
    return data #long, short, PSD

#%%
def CalculateElongEshortPSDWaveforms(waveforms,nSamplesBaseline,nSamplesShort,nSamplesLong):
    baseline = np.mean(waveforms[:,:nSamplesBaseline],axis=1)
    eshort = np.sum(baseline[:,None] - waveforms[:,nSamplesBaseline:nSamplesBaseline+nSamplesShort],axis=1)
    elong = np.sum(baseline[:,None] - waveforms[:,nSamplesBaseline:nSamplesBaseline+nSamplesLong],axis=1)
    PSD = (elong-eshort)/elong
    data = np.hstack((elong.reshape(len(elong),1),eshort.reshape(len(eshort),1),PSD.reshape(len(PSD),1)))
    data = np.hstack((data,waveforms))
    return data #long, short, PSD, waveforms

#%%
def PSDFoM(data):
    n,bins = np.histogram(data[:,2],range=(0,1),bins=1000)

    #Find peaks - initial values used for fit
    peaks, _ = find_peaks(n)
    prominences = peak_prominences(n, peaks)[0]
    maxprominences1 = np.argmax(prominences)
    prominences[maxprominences1] = 0
    maxprominences2 = np.argmax(prominences)

    if (maxprominences1<maxprominences2):
        max1i = peaks[maxprominences1]
        max2i = peaks[maxprominences2]
    else:
        max1i = peaks[maxprominences2]
        max2i = peaks[maxprominences1]

    # Set x values - center of bins
    x = np.empty(n.size)
    for i in range(n.size):
        x[i] = (bins[i]+bins[i+1])/2

    # Define the function to fit to
    def doublegaussian(x, a1, mean1, sigma1, a2, mean2, sigma2,):
        return a1 * np.exp(-((x - mean1)**2 / (2 * sigma1**2))) + a2 * np.exp(-((x - mean2)**2 / (2 * sigma2**2)))

    # Fit data to function
    popt2, _ = curve_fit(doublegaussian, x, n, p0=[n[max1i], bins[max1i], 0.01, n[max2i], bins[max2i], 0.02])

    # Calculate values
    FWHM1 = 2*np.sqrt(2*np.log(2))*popt2[2]
    FWHM2 = 2*np.sqrt(2*np.log(2))*popt2[5]
    FoM = (popt2[4]-popt2[1])/(FWHM1+FWHM2) # FoM = (mean2-mean1)/(FWHM1+FWHM2)
    return FoM

#%%
def PlotPSDandFit(data):
    fig, ax = plt.subplots(constrained_layout=True)
    n,bins,_=plt.hist(data[:,2],range=(0,1),bins=1000,histtype='step')

    #Find peaks - initial values used for fit
    peaks, _ = find_peaks(n)
    prominences = peak_prominences(n, peaks)[0]
    maxprominences1 = np.argmax(prominences)
    prominences[maxprominences1] = 0
    maxprominences2 = np.argmax(prominences)

    if (maxprominences1<maxprominences2):
        max1i = peaks[maxprominences1]
        max2i = peaks[maxprominences2]
    else:
        max1i = peaks[maxprominences2]
        max2i = peaks[maxprominences1]

    # Set x values - center of bins
    x = np.empty(n.size)
    for i in range(n.size):
        x[i] = (bins[i]+bins[i+1])/2

    # Define the function to fit to
    def doublegaussian(x, a1, mean1, sigma1, a2, mean2, sigma2,):
        return a1 * np.exp(-((x - mean1)**2 / (2 * sigma1**2))) + a2 * np.exp(-((x - mean2)**2 / (2 * sigma2**2)))

    # Fit data to function
    popt2, _ = curve_fit(doublegaussian, x, n, p0=[n[max1i], bins[max1i], 0.01, n[max2i], bins[max2i], 0.02])

    # Calculate values
    FWHM1 = 2*np.sqrt(2*np.log(2))*popt2[2]
    FWHM2 = 2*np.sqrt(2*np.log(2))*popt2[5]
    FoM = (popt2[4]-popt2[1])/(FWHM1+FWHM2) # FoM = (mean2-mean1)/(FWHM1+FWHM2)

    # Plot data and fit
    xfit=np.linspace(bins[0],bins[-1],10000)
    plt.plot(xfit,doublegaussian(xfit,*popt2))
    plt.xlabel('PSD')
    plt.ylabel('Counts')
    plt.legend([Line2D([0],[0],color='C0',lw=2),Line2D([0],[0],color='C1',lw=2)],['Data','Fit'])
    plt.xlim(0,0.5)

    textstr = '\n'.join((
    r'$A_1 = %.4f$' % (popt2[0]),
    r'$\mu_1 = %.4f$' % (popt2[1]),
    r'$\sigma_1 = %.4f$' % (popt2[2]),
    r'$A_2 = %.4f$' % (popt2[3]),
    r'$\mu_2 = %.4f$' % (popt2[4]),
    r'$\sigma_2 = %.4f$' % (popt2[5]),
    r'FoM$ = %.4f$' % (FoM)))
    #plt.text(0.75, 0.80, textstr, transform=ax.transAxes, fontsize=10, verticalalignment='top')

#%%
def PlotPSDvsEnergy(data):
    plt.figure(constrained_layout=True)
    plt.hist2d(data[:,0],data[:,2],bins=(500,500),range=([0,np.max(data[:,0])],[0,1]),cmin=1,norm=mpl.colors.LogNorm())
    plt.xlabel('Energy [channels]')
    plt.ylabel('Counts')
    plt.colorbar()

#%% Count number of samples in waveform
def DetermineNumberofSamples(filename):
    file = open(filename, 'r')
    line = file.readlines(1) #headerline
    line = file.readlines(1)
    charcount = 0
    for char in line[0]:
        if char == ';':
            charcount+=1
    file.close()
    nsamples = charcount+1-6
    return nsamples

#%% Plot energy (long) histogram
def PlotEnergy(data,nbins):
    plt.figure(constrained_layout=True)
    plt.hist(data[:,0],bins=nbins,range=(0,np.max(data[:,0])),histtype='step')
    plt.xlabel('Energy (long)')
    plt.ylabel('Counts')

import glob

keywords = ["na22", "na-22", "Na22", "Na-22", "cs137", "Cs137",
            "cs-137", "Cs-137", "AmBe", "ambe"]

#%% Create dependent directories if they do not exist
def dir_check():
    # check if a NPY directory exists already
    npy_exist = os.path.exists("data/DataNPY")
    if npy_exist:
        print("NPY directory does exist")
        NPYfiles = glob.glob("data/DataNPY/*.npy")

        npyfiles = [x.split("\\")[1] for x in NPYfiles]
        npyfiles = [x.split(".npy")[0] for x in npyfiles]

        datasets = [] # na22, cs137, ambe will be the order
        for file in npyfiles:
            datasets.append([x for x in keywords if x in file][0])


    else:
        os.mkdir("data/DataNPY")
        os.mkdir("data/outputs")

        text_datafiles = glob.glob("data/*.txt")
        print(text_datafiles)
        npyfiles = []
        for afile in text_datafiles:
            data = np.loadtxt(afile)
            afile = afile.split('.txt')
            np.save("data/DataNPY/"+afile+'.npy',data)
            npyfiles.append(afile)


        datasets = [] # na22, cs137, ambe will be the order
        for txtfile in text_datafiles:
            datasets.append([x for x in keywords if x in txtfile][0])


    plots_exist = os.path.exists("plots/")
    if plots_exist:
        print("Plot directory does exist")
    else:
        os.mkdir("plots")

    return(datasets,npyfiles)

nakeywords = ["na22", "na-22", "Na22", "Na-22"]
cskeywords = ["cs137", "Cs137", "cs-137", "Cs-137"]
ambekeywords = ["AmBe", "ambe"]

#%%
def automated_routine():
    #%% Fixed values
    nSamplesBaseline = 40
    nSamplesLong = 200
    datasets,filenames = dir_check()

    Na22file = "nofile" ; Cs137file = "nofile" ; AmBefile = "nofile"

    for i in range(len(datasets)):
        if datasets[i] in nakeywords:
            Na22file = filenames[i]
        elif datasets[i] in cskeywords:
            Cs137file = filenames[i]
        elif datasets[i] in ambekeywords:
            AmBefile = filenames[i]

    if Na22file != "nofile":
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Na-22 Source
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Load data
        Na22filepathname = 'data/DataNPY/'+Na22file+'.npy'
        dataNa22 = np.load(Na22filepathname)

        #Plot waveforms (first 10) and baseline, short gate and long gate
        nSamplesShort = 24
        PlotWaveforms(dataNa22[:,3:],10,nSamplesBaseline,nSamplesShort,nSamplesLong)
        plt.title(Na22file)
        plt.savefig("plots/Na22_waveform.png",bbox_inches='tight',dpi=400)
        # Integrate the waveforms to calculate Qlong, Qshort and the PSD
        data = CalculateElongEshortPSD(dataNa22[:,3:],nSamplesBaseline,nSamplesShort,nSamplesLong)
        PlotPSDvsEnergy(data)
        plt.title(Na22file)
        plt.savefig("plots/Na22_psd.png",bbox_inches='tight',dpi=400)
        #Plot energy (long) histogram
        plt.figure(constrained_layout=True)
        n,bins,_=plt.hist(data[:,0],bins=1000,range=(0,36000),histtype='step')
        plt.xlabel('Energy (long) ')
        plt.xlabel('Energy [channels]')
        plt.ylabel('Counts')
        plt.title(Na22file)
        plt.savefig("plots/Na22_spectrum_uncal.png",dpi=400,bbox_inches='tight')
        f = open("data/outputs/"+"Na22_spectrum"+'.txt','w')
        for i in n:
            f.write('{0:.0f}\n'.format(i))
        f.close()
        #%% Find the peaks
        #plt.figure(constrained_layout=True)
        # n,pbins,_ = plt.hist(data[:,0],bins=100,range=(0,36000),histtype='step')
        n,pbins = np.histogram(data[:,0],bins=100,range=(0,36000))
        #plt.xlabel('Energy (long)')
        #plt.ylabel('Counts')
        peaks,prop = find_peaks(n,prominence=300)
        #plt.plot(pbins[peaks],n[peaks],'r*')
        #plt.title(Na22file)
        if len(peaks) < 3:
            raise Exception("Unable to find peaks")

        #%% Analyze 0.511 MeV Compton Edge

        #Initial values
        emin = pbins[peaks[1]]
        emax = pbins[peaks[1]]*1.5
        dmu = 1
        count = 0
        mu_old = 100000

        #Loop unitl mean converges
        #plt.figure(constrained_layout=True)
        while dmu>0.001:
            count+=1 #interation counter
            nbins = int((emax-emin)/20) #determine number of bins
            labelname = 'Iteration: '+str(count)
            #n,bins,_ = plt.hist(data[:,0],bins=nbins,range=(emin,emax),histtype='step',label=labelname)
            n,bins = np.histogram(data[:,0],bins=nbins,range=(emin,emax))
            x = np.empty(n.size)
            for i in range(n.size):
                x[i] = (bins[i]+bins[i+1])/2
            def gaussian(x, a, mean, sigma, c):
                return a * np.exp(-((x - mean)**2 / (2 * sigma**2))) + c
            popt,pcov = curve_fit(gaussian, x, n, p0=[max(n), x[np.argmax(n)], 1000, n[-1]])#, sigma=np.sqrt(n),absolute_sigma=True)
            dmu = abs(mu_old-popt[1])/mu_old

            #Save old values
            mu_old = popt[1]
            emin_old = emin
            emax_old = emax

            #Cacluate new range values
            emin = popt[1]-0.2*popt[2]
            emax = popt[1]+3*popt[2]

            #Break if stuck in loop
            if count==500:
                raise Exception("Unable to converge")
                break

        # Get uncertainty of fit values
        perr = np.sqrt(np.diag(pcov))

        #Calculate R^2
        residuals = n-gaussian(x,*popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((n-np.mean(n))**2)
        r_squared = 1 - (ss_res / ss_tot)

        #Caluclate resolution
        FWHM = 2*np.sqrt(2*np.log(2))*popt[2]
        res = FWHM/popt[1]*100

        #Plot fit
        #xfit=np.linspace(bins[0],bins[-1],10000)
        #plt.plot(xfit,gaussian(xfit,*popt),'r',label='Fit')
        #plt.xlabel('Energy [channels]')
        #plt.ylabel('Counts')
        #plt.legend()
        #plt.legend([Line2D([0],[0],color='C0',lw=2),Line2D([0],[0],color='r',lw=2)],['Data','Fit'])
        #plt.title(Na22file+' 0.511 MeV Compton Edge')

        #Save values
        Na22_511_chmin = emin_old
        Na22_511_chmax = emax_old
        Na22_511_CE = popt[1]+FWHM/2
        Na22_511_CE_unc = np.sqrt((perr[1])**2+(perr[2]*np.sqrt(2*np.log(2)))**2)
        Na22_511_res = res
        Na22_511_r_squared = r_squared

        #%% Analyze 1.062 MeV Compton Edge

        #Initial values
        emin = pbins[peaks[2]]
        emax = pbins[peaks[2]]*1.5
        dmu = 1
        count = 0
        mu_old = 100000

        #Loop unitl mean converges
        #plt.figure(constrained_layout=True)
        while dmu>0.001:
            count+=1 #interation counter
            nbins = int((emax-emin)/40) #determine number of bins
            labelname = 'Iteration: '+str(count)
            #n,bins,_ = plt.hist(data[:,0],bins=nbins,range=(emin,emax),histtype='step',label=labelname)
            n,bins = np.histogram(data[:,0],bins=nbins,range=(emin,emax))
            x = np.empty(n.size)
            for i in range(n.size):
                x[i] = (bins[i]+bins[i+1])/2
            def gaussian(x, a, mean, sigma, c):
                return a * np.exp(-((x - mean)**2 / (2 * sigma**2))) + c
            popt,pcov = curve_fit(gaussian, x, n, p0=[max(n), x[np.argmax(n)], 1000, n[-1]])#, sigma=np.sqrt(n),absolute_sigma=True)
            dmu = abs(mu_old-popt[1])/mu_old

            #Save old values
            mu_old = popt[1]
            emin_old = emin
            emax_old = emax

            #Cacluate new range values
            emin = popt[1]-0.3*popt[2]
            emax = popt[1]+4*popt[2]

            #Break if stuck in loop
            if count==500:
                raise Exception("Unable to converge")
                break

        # Get uncertainty of fit values
        perr = np.sqrt(np.diag(pcov))

        #Calculate R^2
        residuals = n-gaussian(x,*popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((n-np.mean(n))**2)
        r_squared = 1 - (ss_res / ss_tot)

        #Caluclate resolution
        FWHM = 2*np.sqrt(2*np.log(2))*popt[2]
        res = FWHM/popt[1]*100

        #Plot fit
        # xfit=np.linspace(bins[0],bins[-1],10000)
        # plt.plot(xfit,gaussian(xfit,*popt),'r',label='Fit')
        # plt.xlabel('Energy [Channels]')
        # plt.ylabel('Counts')
        # plt.legend()
        # plt.title(Na22file+' 1.275 MeV Compton Edge')

        #Save values
        Na22_1275_chmin = emin_old
        Na22_1275_chmax = emax_old
        Na22_1275_CE = popt[1]+FWHM/2
        Na22_1275_CE_unc = np.sqrt((perr[1])**2+(perr[2]*np.sqrt(2*np.log(2)))**2)
        Na22_1275_res = res
        Na22_1275_r_squared = r_squared

        print("  ---------------------------------")
        print("             Na22 Data")
        print("  ---------------------------------")
        print("  CHMIN, CHMAX:  {:.2f}, {:.2f}\n".format(Na22_511_chmin, Na22_511_chmax))
        print("  511 keV Compton Edge:  {:.2f} {:.2f} {:.2f} {:.2f}\n".format(Na22_511_CE,Na22_511_CE_unc, Na22_511_res, Na22_511_r_squared))
        print("  CHMIN, CHMAX:  {:.2f}, {:.2f}\n".format(Na22_1275_chmin, Na22_1275_chmax))
        print("  1275 keV Compton Edge: {:.2f} {:.2f} {:.2f} {:.2f}".format(Na22_1275_CE, Na22_1275_CE_unc, Na22_1275_res, Na22_1275_r_squared))
        print("  ---------------------------------")

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Cs-137 Source
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Cs137file != "nofile":

        # Load data
        Cs137filepathname = 'data/DataNPY/'+Cs137file+'.npy'
        dataCs137 = np.load(Cs137filepathname)

        #Plot waveforms (first 10) and baseline, short gate and long gate
        nSamplesShort = 24
        PlotWaveforms(dataCs137[:,3:],10,nSamplesBaseline,nSamplesShort,nSamplesLong)
        plt.title(Cs137file)
        plt.savefig("plots/Cs137_waveform.png",dpi=400, bbox_inches='tight')
        # Integrate the waveforms to calculate Qlong, Qshort and the PSD
        data = CalculateElongEshortPSD(dataCs137[:,3:],nSamplesBaseline,nSamplesShort,nSamplesLong)
        PlotPSDvsEnergy(data)
        plt.title(Cs137file)
        plt.savefig("plots/Cs137_psd.png",dpi=400, bbox_inches='tight')
        # Plot energy (long) histogram
        plt.figure(constrained_layout=True)
        n,bins,_=plt.hist(data[:,0],bins=1000,range=(0,36000),histtype='step')
        plt.xlabel('Energy (long)')
        plt.ylabel('Counts')
        plt.title(Cs137file)
        plt.savefig("plots/Cs137_spectrum_uncal.png", dpi = 400, bbox_inches='tight')
        #Output histogram data to file
        f = open("data/outputs/"+"Cs137_spectrum"+'.txt','w')
        for i in n:
            f.write('{0:.0f}\n'.format(i))
        f.close()
        #%% Find the peaks
        # plt.figure(constrained_layout=True)
        #n,pbins,_ = plt.hist(data[:,0],bins=100,range=(0,36000),histtype='step')
        n,pbins = np.histogram(data[:,0],bins=100,range=(0,36000))
        # plt.xlabel('Energy (long)')
        # plt.ylabel('Counts')
        peaks,prop = find_peaks(n,prominence=200)
        # plt.plot(pbins[peaks],n[peaks],'r*')
        # plt.title(Cs137file)

        if len(peaks) < 2:
            raise Exception("Unable to find peaks")

        #%% Analyze 0.662 MeV Compton Edge

        #Initial values
        emin = pbins[peaks[1]]
        emax = pbins[peaks[1]]*1.5
        dmu = 1
        count = 0
        mu_old = 100000

        #Loop unitl mean converges
        #plt.figure(constrained_layout=True)
        while dmu>0.001:
            count+=1 #interation counter
            nbins = int((emax-emin)/20) #determine number of bins
            labelname = 'Iteration: '+str(count)
            #n,bins,_ = plt.hist(data[:,0],bins=nbins,range=(emin,emax),histtype='step',label=labelname)
            n,bins = np.histogram(data[:,0],bins=nbins,range=(emin,emax))
            x = np.empty(n.size)
            for i in range(n.size):
                x[i] = (bins[i]+bins[i+1])/2
            def gaussian(x, a, mean, sigma, c):
                return a * np.exp(-((x - mean)**2 / (2 * sigma**2))) + c
            popt,pcov = curve_fit(gaussian, x, n, p0=[max(n), x[np.argmax(n)], 1000, n[-1]])#, sigma=np.sqrt(n),absolute_sigma=True)
            dmu = abs(mu_old-popt[1])/mu_old

            #Save old values
            mu_old = popt[1]
            emin_old = emin
            emax_old = emax

            #Cacluate new range values
            emin = popt[1]-0.3*popt[2]
            emax = popt[1]+3*popt[2]

            #Break if stuck in loop
            if count==500:
                raise Exception("Unable to converge")
                break

        # Get uncertainty of fit values
        perr = np.sqrt(np.diag(pcov))

        #Calculate R^2
        residuals = n-gaussian(x,*popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((n-np.mean(n))**2)
        r_squared = 1 - (ss_res / ss_tot)

        #Caluclate resolution
        FWHM = 2*np.sqrt(2*np.log(2))*popt[2]
        res = FWHM/popt[1]*100

        #Plot fit
        #xfit=np.linspace(bins[0],bins[-1],10000)
        #plt.plot(xfit,gaussian(xfit,*popt),'r',label='Fit')
        #plt.xlabel('Energy [Channels]')
        #plt.ylabel('Counts')
        #plt.legend()
        #plt.title(Cs137file+' 0.662 MeV Compton Edge')

        #Save values
        Cs137_662_chmin = emin_old
        Cs137_662_chmax = emax_old
        Cs137_662_CE = popt[1]+FWHM/2
        Cs137_662_CE_unc = np.sqrt((perr[1])**2+(perr[2]*np.sqrt(2*np.log(2)))**2)
        Cs137_662_res = res
        Cs137_662_r_squared = r_squared

        #%% Print values to screen
        print("  ---------------------------------")
        print("            Cs137 Data")
        print("  ---------------------------------")
        print("  CHMIN, CHMAX:  {:.2f}, {:.2f}\n".format(Cs137_662_chmin, Cs137_662_chmax))
        print("  Compton Edge:  {:.2f} {:.2f} {:.2f} {:.2f}".format(Cs137_662_CE, Cs137_662_CE_unc, Cs137_662_res, Cs137_662_r_squared))
        print("  ---------------------------------")

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # AmBe Source
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if AmBefile != "nofile":

        # Load data
        AmBefilepathname = 'data/DataNPY/'+AmBefile+'.npy'
        dataAmBe = np.load(AmBefilepathname)

        nSamplesShort = 24
        PlotWaveforms(dataAmBe[:,3:],50,nSamplesBaseline,nSamplesShort,nSamplesLong)
        plt.title(AmBefile)
        plt.savefig("plots/AmBe_waveforms.png", dpi= 400, bbox_inches='tight')
        data = CalculateElongEshortPSDWaveforms(dataAmBe[:,3:],nSamplesBaseline,nSamplesShort,nSamplesLong)

        #%% Energy calibration
        if Na22file != "nofile" and Cs137file != "nofile":
            EMeV = [0.341, 0.477, 1.062]
            Ech = [Na22_511_CE, Cs137_662_CE, Na22_1275_CE]
        if Na22file == "nofile":
            EMeV = [0.477]
            Ech = [Cs137_662_CE]
            print("!!! Warning, light calibration is using only the Cs-137 CE !!!")
            print("!!!   Please provide Na-22 data for a better calibration   !!!")
        if Cs137file == "nofile":
            EMeV = [0.341, 1.062]
            Ech = [Na22_511_CE, Na22_1275_CE]
        # plt.figure(constrained_layout=True)
        # plt.plot(Ech,EMeV,'.')
        # plt.xlabel('Energy [channel]')
        # plt.ylabel('Energy [MeV]')

        pfit = np.polyfit(Ech,EMeV,1)
        #xfit = np.linspace(Ech[0],Ech[-1],100)
        #yfit = np.polyval(pfit,xfit)
        #plt.plot(xfit,yfit)
        #plt.legend(('Data','Fit'))
        #plt.title(dataset)
        # Add calibrated energy to data
        Ecal = np.polyval(pfit,data[:,0])
        data = np.hstack((Ecal.reshape(len(Ecal),1),data))

        #%% Plot energy (calibrated) vs PSD histogram
        plt.figure(constrained_layout=True)
        plt.hist2d(data[:,0],data[:,3],bins=(500,500),range=([0,np.max(data[:,0])],[0,1]),cmin=1,norm=mpl.colors.LogNorm())
        plt.xlabel('Energy [MeVee]')
        plt.ylabel('PSD')
        plt.colorbar()
        plt.title(AmBefile)

        #%% Cut out events between 450 and 550 keV
        datacut = data[data[:,0]>=0.450]
        datacut = datacut[datacut[:,0]<=0.550]

        plt.figure(constrained_layout=True)
        plt.hist2d(data[:,0],data[:,3],bins=(500,500),range=([0,np.max(data[:,0])],[0,1]),cmin=1,norm=mpl.colors.LogNorm())
        plt.xlabel('Energy [MeVee]')
        plt.ylabel('PSD')
        plt.colorbar()
        plt.axvline(0.450,color='tab:red')
        plt.axvline(0.550,color='tab:red')
        plt.ylim(0,0.6)
        plt.title(AmBefile)
        plt.savefig('plots/AmBe_psd.png', bbox_inches='tight',dpi=400)

        #%% Loop through different short gates to determine best PSD FoM
        maxFoM = 0
        maxFoMnSamplesShort = 0
        for nSamplesShort in np.arange(15,30,1):
            dataPSD = CalculateElongEshortPSD(datacut[:,4:],nSamplesBaseline,nSamplesShort,nSamplesLong)
            FoM = PSDFoM(dataPSD)
            if FoM>maxFoM:
                maxFoM = FoM
                maxFoMnSamplesShort = nSamplesShort
        print("            FoM Data")
        print("  ---------------------------------")
        print("  Short Gate, FoM:  {},  {:.3f}".format(maxFoMnSamplesShort,maxFoM))
        print("  ---------------------------------")

        #%% Plot PSD FoM
        nSamplesShort = maxFoMnSamplesShort
        dataPSD = CalculateElongEshortPSD(datacut[:,4:],nSamplesBaseline,nSamplesShort,nSamplesLong)
        PlotPSDandFit(dataPSD)
        plt.title(AmBefile)
        plt.savefig('plots/AmBe_FoM.png',bbox_inches='tight',dpi=400)

        #%% Output data to file
        nSamplesShort = 24
        dataPSD = CalculateElongEshortPSD(data[:,4:],nSamplesBaseline,nSamplesShort,nSamplesLong)
        plt.figure(constrained_layout=True)
        plt.hist2d(dataPSD[:,0],dataPSD[:,2],bins=(500,500),range=([0,np.max(dataPSD[:,0])],[0,1]),cmin=1,norm=mpl.colors.LogNorm())
        #Output histogram data to file
        f = open("data/outputs/"+"AmBe_psd"+'.txt','w')
        f.write('E_long[keVee]\t E_short[keVee]\n')
        for i in range(len(data[:,0])):
            f.write('{:.4f}\t{:.4f}\n'.format(np.polyval(pfit,dataPSD[i,0]),dataPSD[i,2]))
        f.close()

        #d = np.loadtxt("data/outputs/"+AmBefile+'.txt',skiprows=1)
        #plt.figure(constrained_layout=True)
        #d=d[d[:,0]>0]
        #plt.hist2d(d[:,0],d[:,1],bins=(500,500),range=([0,np.max(d[:,0])],[0,1]),cmin=1,norm=mpl.colors.LogNorm())

        #%%
        # cuts = [0,10000,20000,30000,40000,50000,60000,70000,100000]
        # FoMs = []

        # plt.figure(constrained_layout=True)
        # for i in range(1,len(cuts)-1):
        #     datacut = data[data[:,0]>cuts[i]]
        #     datacut = datacut[datacut[:,0]<cuts[i+1]]
        #     FoM = PSDFoM(datacut)
        #     print(FoM)
        #     plt.errorbar((cuts[i]+cuts[i+1])/2,FoM,xerr=(cuts[i]+cuts[i+1])/2-cuts[i],marker='.',color='C0')

        #%% Plot waveforms
        # plt.rcParams.update({'font.size': 12})

        # plt.figure(constrained_layout=True)
        # plt.hist2d(data[:,0],data[:,3],bins=(500,500),range=([0,np.max(data[:,0])],[0,1]),cmin=1,norm=mpl.colors.LogNorm())
        # plt.xlabel('Energy [MeVee]')
        # plt.ylabel('PSD')
        # plt.colorbar()
        # plt.title(AmBefile)

        # dataecut = data[data[:,0]>=0.9]
        # dataecut = dataecut[dataecut[:,0]<=1.1]

        # gdata = dataecut[dataecut[:,3]>=0.1]
        # gdata = gdata[gdata[:,3]<=0.14]

        # ndata = dataecut[dataecut[:,3]>=0.19]
        # ndata = ndata[ndata[:,3]<=0.23]

        # navewave = np.mean(ndata[:,4:],axis=0)
        # gavewave = np.mean(gdata[:,4:],axis=0)

        # time = np.arange(0,496,2)
        # plt.figure(constrained_layout=True)
        # plt.plot(time,gavewave,linewidth=3,label='Gamma')
        # plt.plot(time,navewave,linewidth=3,label='Neutron')
        # plt.ylabel('ADC [channel]')
        # plt.xlabel('Time [ns]')
        # plt.legend()
        # plt.xlim(75,250)
        # plt.savefig('AFIT101waveforms.pdf')

automated_routine()
