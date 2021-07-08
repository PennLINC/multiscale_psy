%%% Merge all subject/scale-level FC matrices into one struct for further proc

% to become 2 to 30 when stuff finishes running someday
Krange=2:25;

subjs=load('/cbica/projects/pinesParcels/data/bblids.txt');
psySubjs=load('/cbica/projects/pinesParcels/data_psy/NewPsySubjs.txt');

% initialize structures to be filled
for i=2:max(Krange)
	ind_mats{i}=zeros(i,i,(length(subjs)+length(psySubjs)));
	gro_mats{i}=zeros(i,i,(length(subjs)+length(psySubjs)));
	bTS_indmats{i}=zeros(i,i,(length(subjs)+length(psySubjs)));
end

% fill in with each subj, note K=1 is empty for each
for s=1:length(subjs)
	% load in data
	fp=['/cbica/projects/pinesParcels/data/CombinedData/' num2str(subjs(s)) '/fc_metrics.mat']
	fcmets=load(fp);	
	% throw it in 
	for k=2:max(Krange)
		ind_mats{k}(:,:,s)=fcmets.subjmats(k).Khouse;
		gro_mats{k}(:,:,s)=fcmets.subjmats(k).GKhouse;
		bTS_indmats{k}(:,:,s)=fcmets.subjmats(k).K_bTS_house;
	end
end

% stack new psy subjs on top
for s=1:length(psySubjs)
        % load in data
        fp=['/cbica/projects/pinesParcels/data_psy/CombinedData/' num2str(psySubjs(s)) '/fc_metrics.mat']
        fcmets=load(fp);
        % throw it in 
        for k=2:max(Krange)
                ind_mats{k}(:,:,s+length(subjs))=fcmets.subjmats(k).Khouse;
                gro_mats{k}(:,:,s+length(subjs))=fcmets.subjmats(k).GKhouse;
                bTS_indmats{k}(:,:,s+length(subjs))=fcmets.subjmats(k).K_bTS_house;
        end
end

% save the cell struct array frakenmatrices
save('/cbica/projects/pinesParcels/results_psy/aggregated_data/ind_conmats_allscales_allsubjs.mat','ind_mats');
save('/cbica/projects/pinesParcels/results_psy/aggregated_data/gro_conmats_allscales_allsubjs.mat','gro_mats');
save('/cbica/projects/pinesParcels/results_psy/aggregated_data/bts_conmats_allscales_allsubjs.mat','bTS_indmats');
