import mat73
import numpy as np
from scipy.signal import hilbert
from scipy import signal
from scipy.signal import find_peaks
from scipy.signal import butter

# initialize dataframe for subject ID, bottom up and top down instances, length of wave in TRs
dataframe=np.zeros([790,5])

# load in subject data
subjs=open('/cbica/projects/pinesParcels/data/bblids.txt')
subjs=subjs.read()
subjs=subjs.splitlines()
psySubjs=open('/cbica/projects/pinesParcels/data_psy/NewPsySubjs.txt')
psySubjs=psySubjs.read()
psySubjs=psySubjs.splitlines()

# set info on sampling frequency (TR) and time series length
fs=.33
# nyquist limit
Nyq=fs/2
# adopt desired bandpass to nyquist limit/butterworth filtering
desired_low=.01/Nyq
desired_high=.05/Nyq
# bandpass filter to be used
b, a = butter(2,[desired_low,desired_high],'band')

for s in range(len(subjs)):
	# load data
	subjFP='/cbica/projects/pinesParcels/data/SingleParcellation/SingleParcel_1by1_kequal_4/Sub_' + subjs[s] + '/IndividualParcel_Final_sbj1_comp4_alphaS21_1_alphaL10_vxInfo1_ard0_eta0/final_UV.mat'
	subjData=mat73.loadmat(subjFP,'r')
	# extract basis time series
	sU=np.array(subjData['U'])
	# arrange in order of hierarchy
	TMOrder_4=[1,0,2,3]
	sU_TMOrdered=sU[0,0:554,TMOrder_4]
	# apply bandpass filter
	for x in range(4):
		ts=signal.filtfilt(b,a,sU_TMOrdered[x,:])
		# adding time series mean to take out negativity
		sU_TMOrdered[x,:]=ts+np.mean(sU_TMOrdered[x,:])	
	
	# extract motor peaks
	x=0
	# transform into phase of time series
	z=hilbert(sU_TMOrdered[x,0:554])
	ts=np.unwrap(np.angle(z))
	# find peaks in phase, at least 4 TRs apart
	Motor_peaks, _ = find_peaks(ts, distance=4)
	# threshold out peaks with negative phase value
	Motor_peaks = Motor_peaks[ts[Motor_peaks]>0]	
	# extract FP peaks
	x=2
	# transform into phase of time series
	z=hilbert(sU_TMOrdered[x,0:554])
	ts=np.unwrap(np.angle(z))
	# find peaks in phase, at least 4 TRs apart
	FP_peaks, _ = find_peaks(ts, distance=4)
	# threshold out peaks with negative phase value
	FP_peaks = FP_peaks[ts[FP_peaks]>0]
	# extract DM peaks
	x=3
	# transform into phase of time series
	z=hilbert(sU_TMOrdered[x,0:554])
	ts=np.unwrap(np.angle(z))
	# find peaks in phase, at least 3 TRs apart
	DM_peaks, _ = find_peaks(ts, distance=4)
	# threshold out peaks with negative phase value
	DM_peaks = DM_peaks[ts[DM_peaks]>0]
	#Bottom up prop count
	BU=0
	# wave duration count
	BUWD=0
	# for each motor peak
	for P in range(len(Motor_peaks)):
    		# get timepoint of motor peak
		tp=Motor_peaks[P]
    		# are there DM and FP peaks past this motor peak?
		if tp < max(DM_peaks) and tp < max(FP_peaks):
        		### find TP of next FP peak
        		# all FP peaks > mot peak
			AllFP_G=FP_peaks[FP_peaks>tp]
        		# immediately proceeding peak
			NextFP_P=AllFP_G[0]
        		### find TP of next DM peak
			AllDM_G=DM_peaks[DM_peaks>tp]
        		# immediately proceeding peak
			NextDM_P=AllDM_G[0]    
        		# does it pass the sniff test? (DM TP > FP TP, full procession < 17 TRs)
			if ((NextDM_P > NextFP_P) and ((NextDM_P - tp) < 17)):
				# add instance to iterative count of bottom-up props
				BU=BU+1
				# add duration of wave, change to units of seconds from TRs
				BUWD=BUWD+((NextDM_P-tp)*3)
	
	# top down prop count
	TD=0
	# wave duration count
	TDWD=0
	# for each DM peak
	for P in range(len(DM_peaks)):
    		# get timepoint of DM peak
		tp=DM_peaks[P]
    		# are there DM and FP peaks past this motor peak?
		if tp < max(Motor_peaks) and tp < max(FP_peaks):
        		### find TP of next FP peak
        		# all FP peaks > mot peak
			AllFP_G=FP_peaks[FP_peaks>tp]
        		# immediately proceeding peak
			NextFP_P=AllFP_G[0]
        		### find TP of next DM peak
			AllMO_G=Motor_peaks[Motor_peaks>tp]
        		# immediately proceeding peak
			NextMO_P=AllMO_G[0]    
        		# does it pass the sniff test? (DM TP > FP TP, full procession < 17 TRs)
			if ((NextMO_P > NextFP_P) and ((NextMO_P - tp) < 17)):
            			# add instance to iterative count of bottom-up props
				TD=TD+1
				# add duration of wave, change to units of seconds from TRs
				TDWD=TDWD+((NextMO_P-tp)*3)
	
	# add subject and instances of BD and TU to dataframe
	dataframe[s,0]=subjs[s]
	dataframe[s,1]=BU
	dataframe[s,2]=TD
	# get average wave duration for BU and TD
	a_BUWD=BUWD/BU
	dataframe[s,3]=a_BUWD
	a_TDWD=TDWD/TD
	dataframe[s,4]=a_TDWD

# for psych. subjs
for s in range(len(psySubjs)):
	# load data
	subjFP='/cbica/projects/pinesParcels/data_psy/SingleParcellation/SingleParcel_1by1_kequal_4/Sub_' + psySubjs[s] + '/IndividualParcel_Final_sbj1_comp4_alphaS21_1_alphaL10_vxInfo1_ard0_eta0/final_UV.mat'
	subjData=mat73.loadmat(subjFP,'r')
	# extract basis time series
	sU=np.array(subjData['U'])
	# arrange in order of hierarchy
	TMOrder_4=[1,0,2,3]
	sU_TMOrdered=sU[0,0:554,TMOrder_4]
        # apply bandpass filter
	for x in range(4):
		ts=signal.filtfilt(b,a,sU_TMOrdered[x,:])
		# adding time series mean to take out negativity
		sU_TMOrdered[x,:]=ts+np.mean(sU_TMOrdered[x,:])

	# extract motor peaks
	x=0
	# transform into phase of time series
	z=hilbert(sU_TMOrdered[x,0:554])
	ts=np.unwrap(np.angle(z))
	# find peaks in phase, at least 4 TRs apart
	Motor_peaks, _ = find_peaks(ts, distance=4)
	# threshold out peaks with negative phase value
	Motor_peaks = Motor_peaks[ts[Motor_peaks]>0]	
	# extract FP peaks
	x=2
	# transform into phase of time series
	z=hilbert(sU_TMOrdered[x,0:554])
	ts=np.unwrap(np.angle(z))
	# find peaks in phase, at least 4 TRs apart
	FP_peaks, _ = find_peaks(ts, distance=4)
	# threshold out peaks with negative phase value
	FP_peaks = FP_peaks[ts[FP_peaks]>0]
	# extract DM peaks
	x=3
	# transform into phase of time series
	z=hilbert(sU_TMOrdered[x,0:554])
	ts=np.unwrap(np.angle(z))
	# find peaks in phase, at least 3 TRs apart
	DM_peaks, _ = find_peaks(ts, distance=4)
	# threshold out peaks with negative phase value
	DM_peaks = DM_peaks[ts[DM_peaks]>0]
	#Bottom up prop count
	BU=0
        # wave duration count
	BUWD=0
	# for each motor peak
	for P in range(len(Motor_peaks)):
		# get timepoint of motor peak
		tp=Motor_peaks[P]
    		# are there DM and FP peaks past this motor peak?
		if tp < max(DM_peaks) and tp < max(FP_peaks):
        		### find TP of next FP peak
        		# all FP peaks > mot peak
			AllFP_G=FP_peaks[FP_peaks>tp]
        		# immediately proceeding peak
			NextFP_P=AllFP_G[0]
        		### find TP of next DM peak
			AllDM_G=DM_peaks[DM_peaks>tp]
        		# immediately proceeding peak
			NextDM_P=AllDM_G[0]    
        		# does it pass the sniff test? (DM TP > FP TP, full procession < 17 TRs)
			if ((NextDM_P > NextFP_P) and ((NextDM_P - tp) < 17)):
            			# add instance to iterative count of bottom-up props
				BU=BU+1
                                # add duration of wave, change to units of seconds from TRs
				BUWD=BUWD+((NextDM_P-tp)*3)	
	
	#top down prop count
	TD=0
        # wave duration count
	TDWD=0
	# for each DM peak
	for P in range(len(DM_peaks)):
    		# get timepoint of DM peak
		tp=DM_peaks[P]
    		# are there DM and FP peaks past this motor peak?
		if tp < max(Motor_peaks) and tp < max(FP_peaks):
        		### find TP of next FP peak
        		# all FP peaks > mot peak
			AllFP_G=FP_peaks[FP_peaks>tp]
        		# immediately proceeding peak
			NextFP_P=AllFP_G[0]
        		### find TP of next DM peak
			AllMO_G=Motor_peaks[Motor_peaks>tp]
        		# immediately proceeding peak
			NextMO_P=AllMO_G[0]    
        		# does it pass the sniff test? (DM TP > FP TP, full procession < 17 TRs)
			if ((NextMO_P > NextFP_P) and ((NextMO_P - tp) < 17)):
            			# add instance to iterative count of bottom-up props
				TD=TD+1
                                # add duration of wave, change to units of seconds from TRs
				TDWD=TDWD+((NextMO_P-tp)*3)
	
	# add subject and instances of BD and TU to dataframe
	dataframe[s+692,0]=psySubjs[s]
	dataframe[s+692,1]=BU
	dataframe[s+692,2]=TD
        # get average wave duration for BU and TD
	a_BUWD=BUWD/BU
	dataframe[s+692,3]=a_BUWD
	a_TDWD=TDWD/TD
	dataframe[s+692,4]=a_TDWD

np.savetxt('BUTD.csv',dataframe,delimiter=",")
