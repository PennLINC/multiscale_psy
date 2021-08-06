# Normalization, binning of vertexwise time series, wave property saveout 
import scipy
import nibabel as nb
import numpy as np
from scipy import stats
from scipy import signal
from scipy.signal import find_peaks
from scipy.signal import butter
from scipy import signal
from numpy import genfromtxt
import sys
import sklearn
from sklearn import linear_model

#### load in nonSubj TS data
# load in principal gradient
PG=nb.load('/cbica/projects/abcdfnets/data/hcp.gradients.dscalar.nii')

# Subject is set to the passed argument
subj = sys.argv[1]

# all the scan types
tasks=['rest','SST','nback','mid'];

# parent filepath
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'
# child (output) filepath
childfp='/cbica/projects/abcdfnets/results/' + str(subj) + '/'
# for each task
for T in range(len(tasks)):
	# initialize big array for distribution of all PG delay correlations
	PGD_arr=[]
	# load in time series (masked, bp filtered)
	filepath=parentfp + str(subj) + '_p2mm_masked_filtered_' + tasks[T] + '.dtseries.nii'
	subjData=nb.load(filepath)
	# load in time series
	procTS=subjData.dataobj
	# normalize to mean and SD
	Avg=np.mean(procTS,axis=0)
	SD=np.std(procTS,axis=0)
	procTS=(procTS-Avg)/SD
	# convert TS to numpy array to allow for indexing
	procTS=np.array(procTS)
	# initialize empty array for gradient bins
	procTS_bins=np.zeros((len(procTS),100))
	# bin vertices at same position on gradient
	for b in range(100):
		gradPrctile=np.percentile(PG.dataobj[0,:],b)
		gradPrctile_upper=np.percentile(PG.dataobj[0,:],b+1)
		# index of vertices belonging to this percentile
		boolean_of_interest=np.logical_and(PG.dataobj[0,:] > gradPrctile, PG.dataobj[0,:] < gradPrctile_upper)
		PGindices=np.nonzero(boolean_of_interest)
		PGindices_array=np.array(PGindices)
		# initialize array of all these vertices to average over
		initPGbin=procTS[:,PGindices_array]
		# average signal over this bin 
		meanSig=np.mean(initPGbin,axis=2)
		meanSig=meanSig[:,0]
		# plop into ProcTS_bins
		procTS_bins[:,b]=meanSig
	
	# load in global signal
	GSfFP=parentfp + str(subj) + '_p2mm_masked_filtered_' + tasks[T] + '_GS.csv'
	GS=np.genfromtxt(GSfFP,delimiter=",")
	# normalize GS
	GAvg=np.mean(GS)
	GSD=np.std(GS)
	GS=((GS-GAvg)/GSD)
	# calculate GS troughs with negative find_peaks
	GS_troughs, _ = find_peaks(-GS, distance=8)
	# make an array for each percentile bin and each trough
	troughsNum=len(GS_troughs)-1
	delayMatrix=np.zeros((100,troughsNum))
	# and a relative magnitude matrix
	magMatrix=np.zeros((100,troughsNum))
	# for each trough-trough interval, find peak of bin timeseries
	for t in range(troughsNum):
		tstart=GS_troughs[t]
		tend=GS_troughs[t+1]
		# get GS peak here
		GS_peak, _ = find_peaks(GS[tstart:tend],distance=(tend-tstart))
		for b in range(100):
			# isolate time series sequence
			iso_ts=procTS_bins[tstart:tend,b]
			# find peak in this bin (set min distance to be temporal width of bin)
			peak, _ =find_peaks(iso_ts,distance=(tend-tstart))
			# determine distance from GS peak
			distanceFGSP=peak-GS_peak
			# if peak exists, add to matrix
			if ((len(peak) !=0) and (len(GS_peak) !=0)):
				delayMatrix[b,t]=distanceFGSP
				# record magnitude of normalized signal as point of peak
				magMatrix[b,t]=iso_ts[peak]
			else:
				delayMatrix[b,t]=999
				magMatrix[b,t]=999
	
	# extract columns where each percentile bin had a peak (MORE EXCLUSIVE THAN CERBCORT PAPER)
	delayMatrix_thresh = delayMatrix[:,np.all(delayMatrix != 999, axis=0)]
	# calculate distribution of correlations of (PG location, trough offset) for each interval
	CorDistr=np.zeros((1,delayMatrix_thresh.shape[1]))
	# calculate relative magnitude of wave over its course in units of normalized signal
	Wslopes=np.zeros((1,delayMatrix_thresh.shape[1]))
	for i in range(delayMatrix_thresh.shape[1]):
		CorDistr[0,i], _ =stats.pearsonr(delayMatrix_thresh[:,i],np.arange(100))
		# set x to relative magnitude
		my_x = magMatrix[:,i]
		# set y to position in gradient
		my_y = np.arange(100)
		slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(my_x, my_y)
		Wslopes[0,i] = slope

	# now that we have PG*Delay correlation, delineate the speed in TRs of each wave
	Wspeeds=np.zeros((1,delayMatrix_thresh.shape[1]))
	for i in range(delayMatrix_thresh.shape[1]):
		Wspeeds[0,i] = max(delayMatrix_thresh[:,i])-min(delayMatrix_thresh[:,i])
	
	# now that we have the speed in TRs of each wave, delineate the slope in relative magnitude over TRs
	Worigin=np.zeros((1,delayMatrix_thresh.shape[1]))
	for i in range(delayMatrix_thresh.shape[1]):
		# origin estimation
		val, idx = min((val, idx) for (idx, val) in enumerate(delayMatrix_thresh[:,i]))
		Worigin[0,i]=idx
		# get linear slope of relative magnitude over wave duration
				
	# save out GW # x 4 matrix for this subj for this task: PG*delay corr, Speed, slope, and origin for each GW
	saveoutMat=np.zeros((4,delayMatrix_thresh.shape[1]))
	saveoutMat[0,:]=CorDistr
	saveoutMat[1,:]=Wspeeds
	saveoutMat[2,:]=Worigin
	saveoutMat[3,:]=Wslopes
	saveFN=childfp + str(subj) + '_' + str(tasks[T]) + '_waveProps.csv'
	np.savetxt(saveFN,saveoutMat,delimiter=",")

