# Temporal filtering, normalization, and binning of vertexwise time series

import nibabel as nb
import numpy as np
from scipy.signal import hilbert, chirp
from scipy import stats
from scipy import signal
from scipy.signal import find_peaks
from scipy.signal import butter
from scipy import signal
import pandas as pd
from numpy import genfromtxt

#### construct bandpass filter, load in nonSubj TS data
# sampling frequency
fs=1.389
# nyquist limit
Nyq=fs/2
# adopt desired bandpass to nyquist limit/butterworth filtering
desired_low=.01/Nyq
desired_high=.1/Nyq
# butterworth filter
b, a = butter(2,[desired_low,desired_high],'band')
# load in principal gradient
PG=nb.load('/cbica/projects/pinesParcels/GS/hcp.gradients.dscalar.nii')
# load in subject list
subjs=open('/cbica/projects/pinesParcels/GS/fn.txt')
subjs=subjs.read()
subjs=subjs.splitlines()
# initialize big array for distribution of all PG delay correlations
PGD_arr=[]

#### For each subject
### (for a few subjects first)
for s in range(5):
	print(subjs[s])
	# load in time series
	filepath='/cbica/home/bertolem/xcp_hcp/fmriprepdir/' + str(subjs[s]) + '/func/' + str(subjs[s]) + '_task-REST1_acq-RL_space-fsLR_den-91k_bold.dtseries.nii'
	subjData=nb.load(filepath)
	# initialize empty array to input proc'ed data
	ProcTS=np.zeros((1200,91282))
	# for each vertex
	for x in range(91282):
		# bandpass filter subject data
		ProcTS[:,x]=signal.filtfilt(b,a,subjData.dataobj[:,x])
		# normalize to mean and SD
		Avg=np.mean(ProcTS[:,x])
		SD=np.std(ProcTS[:,x])
		ProcTS[:,x]=((ProcTS[:,x]-Avg)/SD)
		# end for each vertex
	
	print('done with vertex filt/norm')
	# initialize empty array to input proc'ed bin'ed data
	ProcTS_bins=np.zeros((1200,100))
	# bin vertices at same position on gradient
	for b in range(100):
		gradPrctile=np.percentile(PG.dataobj[0,:],b)
		gradPrctile_upper=np.percentile(PG.dataobj[0,:],b+1)
		# index of vertices belonging to this percentile
		boolean_of_interest=np.logical_and(PG.dataobj[0,:] > gradPrctile, PG.dataobj[0,:] < gradPrctile_upper)
		PGindices=np.nonzero(boolean_of_interest)
		PGindices_array=np.array(PGindices)
		# initialize array of all these vertices to average over
		initPGbin=ProcTS[:,PGindices_array]
		# average signal over this bin 
		meanSig=np.mean(initPGbin,axis=2)
		meanSig=meanSig[:,0]
		# plop into ProcTS_bins
		ProcTS_bins[:,b]=meanSig
		# end for each bin

	print('done with PG binning, signal averaging')
	# load in global signal
	ConfFP='/cbica/home/bertolem/xcp_hcp/fmriprepdir/' + str(subjs[s]) + '/' + str(subjs[s]) + '_task-REST2_acq-RL_desc-confounds_timeseries.tsv'
	Conf=np.genfromtxt(ConfFP,delimiter="\t",names=True)
	GS=Conf['global_signal']
	# filter the global signal (bandpass)
	GS_filt=np.zeros(GS.shape)
	GS_filt=signal.filtfilt(b,a,GS)
	# normalize GS
	GAvg=np.mean(GS_filt)
	GSD=np.std(GS_filt)
	GS_filt=((GS_filt-GAvg)/GSD)
	# calculate GS troughs with negative find_peaks
	GS_filt_troughs, _ = find_peaks(-GS_filt, distance=8)
	# make an array for each percentile bin and each trough
	troughsNum=len(GS_filt_troughs)-1
	delayMatrix=np.zeros((100,troughsNum))
	# for each trough-trough interval, find peak of bin timeseries
	for t in range(troughsNum):
		tstart=GS_filt_troughs[t]
		tend=GS_filt_troughs[t+1]
		# get GS peak here
		GS_peak, _ = find_peaks(GS_filt[tstart:tend],distance=(tend-tstart))
		for b in range(100):
			# isolate time series sequence
			iso_ts=ProcTS_bins[tstart:tend,b]
			# find peak in this bin (set min distance to be temporal width of bin)
			peak, _ =find_peaks(iso_ts,distance=(tend-tstart))
			# determine distance from GS peak
			distanceFGSP=peak-GS_peak
			# if peak exists, add to matrix
			if len(peak) !=0:
				delayMatrix[b,t]=distanceFGSP
			else:
				delayMatrix[b,t]=999

	print('done with delay matrix construction')
	# extract columns where each percentile bin had a peak
	delayMatrix_thresh = delayMatrix[:,np.all(delayMatrix != 999, axis=0)]
	# calculate distribution of correlations of (PG location, trough offset) for each interval
	CorDistr=np.zeros((1,delayMatrix_thresh.shape[1]))
	for i in range(delayMatrix_thresh.shape[1]):
		CorDistr[0,i], _ =stats.pearsonr(delayMatrix_thresh[:,i],np.arange(100))

	PGD_arr.append(CorDistr)
	print(subjs[s])

np.savetxt('PGD_arr.csv',PGD_arr,delimiter=",")
