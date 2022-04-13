import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from gaussian import gauss_exp_fit
from scipy.signal import find_peaks, peak_prominences


sns.set_theme()
sns.set_style('ticks')
##------------------------------------------------------------------------
## Step 1. Read in the gamma data from the Am241 source
##------------------------------------------------------------------------
path = '2022_04_04_AFITB10JJM_Am241/UNFILTERED/'
df = pd.read_csv(path+"am241_spectrum.csv" , sep = ';')
##------------------------------------------------------------------------
## Step 2. Bin the spectrum to generate the CoMPASS equivalent spectrum
##------------------------------------------------------------------------
# 2a. bin the data, reduce plays around with the binning size
reduce = 8
nbins = int(4096/reduce) # playing around with the binning

gdata,gbins = np.histogram(df["ENERGY"], nbins)#; gbins = gbins[:-1] # bin the data
#-------------------------------------------------------------------------
# Step 3. Locate the Am241 photopeak and apply a best-fit to it
#-------------------------------------------------------------------------
peaks, _ = find_peaks(gdata) # locates the peaks in the binned spectrum
prominences = peak_prominences(gdata, peaks)[0]
maxprominences1 = np.argmax(prominences) # pulls out the index of the 1st peak
prominences[maxprominences1] = 0
maxprominences2 = np.argmax(prominences) # pulls out the index of the 2nd peak

centroid = peaks[maxprominences2] # 1st peak is low channel noise so we want 2nd

# here we isolate the Am241 photopeak
cut = int(centroid - (reduce*2)) ; stop = int(centroid + (reduce*3))

params = [gdata[centroid], gbins[centroid], 100, 100, 0] # initial guesses
g,u,s = gauss_exp_fit(gbins[cut:stop],gdata[cut:stop],params)
print("Centroid Value: {:.3f}".format(u))
#-------------------------------------------------------------------------
# Step 3. Plot the full Am241 spectrum
#-------------------------------------------------------------------------
fig,ax=plt.subplots(figsize=(4,3))
plt.hist(df["ENERGY"],nbins,histtype='step')
plt.ylabel("Counts") ; plt.xlabel("Light (channels)")

ax.set_xticks([1000,3000,5000,7000],minor=True)
ax.set_yticks([2500,7500,12500,17500],minor=True)

ax.tick_params(axis='x', direction='in',top=True)
ax.tick_params(axis='x', which='minor',direction='in',top=True)
ax.tick_params(axis='y', direction='in',right=True)
ax.tick_params(axis='x', which='minor',direction='in',right=True)
ax.tick_params(axis='y', which='minor',direction='in',right=True)

plt.xlim(0,1500)


plt.savefig("plots/uncal.png",dpi=300,bbox_inches='tight')
plt.close()
##------------------------------------------------------------------------
## Step 4. Plot the Am241 photopeak with the Gaussian + exp bkgrnd fit
##------------------------------------------------------------------------
fig,ax=plt.subplots(figsize=(4,3))
plt.hist(df["ENERGY"],nbins,histtype='step')
plt.xlim(150,650)
plt.ylim(0,1000)

# x is the bin centers
x = np.array([(gbins[i]+gbins[i+1])/2 for i in range(len(gbins)-1)])
plt.plot(x[cut:stop],g,color='r',linewidth=0.85)

ax.set_xticks([150,200,250,350,450,550,650],minor=True)
ax.set_yticks([100,300,500,700,900],minor=True)

ax.tick_params(axis='x', direction='in',top=True)
ax.tick_params(axis='x', which='minor',direction='in',top=True)
ax.tick_params(axis='y', direction='in',right=True)
ax.tick_params(axis='x', which='minor',direction='in',right=True)
ax.tick_params(axis='y', which='minor',direction='in',right=True)

#plt.legend(["Eq. 2 Best Fit","Signal"])

plt.ylabel("Counts") ; plt.xlabel("Energy (channels)") ; plt.tight_layout()
plt.savefig("plots/calibrated_gamma_spectrum.png",dpi=300,bbox_inches='tight')
