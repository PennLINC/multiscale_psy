% load in subj list
subjs=load('/cbica/projects/pinesParcels/data_psy/PsychSubjsScanId.txt');
% bblid in scanid order for matching to time series
bblids=load('/cbica/projects/pinesParcels/data_psy/PsychSubjsBblid_inScanIdOrder.txt');
% initializ bigmat
bigmat=zeros(length(subjs),492,14);
% colnames
colNames=cell(492,14);
% subcort list
subcort_list=['R_Accumbens_Area', 'L_Accumbens_Area', 'R_Amygdala', 'L_Amygdala', 'R_Caudate', 'L_Caudate', 'R_Hippocampus', 'L_Hippocampus', 'R_Pallidum', 'L_Pallidum', 'R_Putamen', 'L_Putamen', ' R_Thalamus_Proper', 'L_Thalamus_Proper'];
% populate with each subject except for 707, who is missing
%for i=setdiff(1:length(subjs),707);
for i=1:2
	% load in subcort TS
	subCortTSfp=['/cbica/projects/pinesParcels/data_psy/Subcort/' num2str(subjs(i)) 'rs_nb_ts.csv'];
	subCortTS=load(subCortTSfp);
	% initialize subcort-cort connectivity vector
	for K=2:30
		% get bblid
		bblid=bblids(i);
		% load in basis TS
		K_Folder = ['/cbica/projects/pinesParcels/data/SingleParcellation/SingleParcel_1by1_kequal_' num2str(K) '/Sub_' num2str(bblids(i))];
		K_part_subj =[K_Folder '/IndividualParcel_Final_sbj1_comp' num2str(K) '_alphaS21_1_alphaL10_vxInfo1_ard0_eta0/final_UV.mat'];
		% check if it exists in 693 dir, look in 790 if not
		if ~exist(K_part_subj)
		K_Folder = ['/cbica/projects/pinesParcels/data_psy/SingleParcellation/SingleParcel_1by1_kequal_' num2str(K) '/Sub_' num2str(bblids(i))]; 
                K_part_subj =[K_Folder '/IndividualParcel_Final_sbj1_comp' num2str(K) '_alphaS21_1_alphaL10_vxInfo1_ard0_eta0/final_UV.mat'];
		end
		% load in partition	
		subj_part=load(K_part_subj);
		% obtain U
		subj_cortTS=subj_part.U{:};
		% truncate cort TS
		subj_cortTS_trunc=subj_cortTS(1:351,:);
		% remove 6 volumes of nback (ASSUMPTION)
		subj_cortTS_trunc(121:126,:)=[];
		% correlation of subcorts and cort
		for N=1:K
			NetCors_w_SubC=corr(subCortTS,subj_cortTS_trunc(:,N));
			bigmat(i,((K+N)-2),:)=NetCors_w_SubC;
			% col-naming as proximate to feature recording as possible, but redundant
			for S=1:14
				colName=['ind_bw_FC_scale' num2str(K) '_Net' num2str(N) '_Sub' subcort_list(S)];
				colNames(((K+N)-2),S)=cellstr(colName);
			end
		end
	end	
	% correlation of subcorts
	for s=1:14
		bigmat(i,465:478,s)=corr(subCortTS(:,s),subCortTS);
		%colnames
		for ss=1:14
			colName=['ind_bw_FC_rsnb_' subcort_list(s) '_' subcort_list(ss)];
			colNames((465+ss),s)=cellstr(colName);
		end
		% load in emo subcort cor
		emo_subCortFCfp=['/cbica/projects/pinesParcels/data_psy/Subcort/' num2str(subjs(i)) 'emoID_fc.csv'];
        	emo_subCortFC=load(emo_subCortFCfp);
		% place row of correlation matrix as emoSub-emoSub
		bigmat(i,479:492,s)=emo_subCortFC(s,:);
	        % colnames
		for ss=1:14
                        colName=['ind_bw_FC_emo_' subcort_list(s) '_' subcort_list(ss)];
                        colNames((479+ss),s)=cellstr(colName);
                end
	end
% print subject
i
end
%%% unwrap
% initialize output csv,+1 for colnames
df=cell(length(subjs)+1,(14*492));
% populate as ascending subcort
for s=1:14
	% each subcort ROI will have 492 features: start each subcortical feature vector 492 * (n - 1), end 492 later  
	% 2:length(subjs)+1 so we can save the first row for column names
	df(2:(length(subjs)+1),((492*(s)-491)):(492*(s)))=num2cell(bigmat(1:(length(subjs)),1:492,s));
	% column names
	% each scale
	for K=2:30
		% each network
		for N=1:K
			% set in column name
			% row one for whole df, start each subcortical colName vector at 492 * (n - 1)
			df(1,((492*(s)-491)+((K+N)-3)):(492*(s)))=colNames((K+N)-3,s);
		end
	end
% name subcortical-subcortical columns
df(1,(((492*(s)-491))+((K+N)-2)):(492*(s)))=colNames(465:492,s);
% this iter's colnames
This_Iterations_ColNames=df(1,((492*(s-1)):(492*(s1))))
end	
% save
writetable(cell2table(df),'~/results_psy_master_fcfeats_Sub.csv')
