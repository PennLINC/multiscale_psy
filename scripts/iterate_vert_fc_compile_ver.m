%iterate over subjects to get within and between network connectivities for each scale (K)

% add needed paths
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
ProjectFolder = '/cbica/projects/pinesParcels/data/SingleParcellation';
ProjectFolder_psy = '/cbica/projects/pinesParcels/data_psy/SingleParcellation';
% What is K range to iterate over?
% I'll tell you hwhat
Krange=2:30;
% Read in subjects list
subjs=load('/cbica/projects/pinesParcels/data_psy/NewPsySubjs.txt');
% read in group partitions
group_parts=load([ProjectFolder '/SingleAtlas_Analysis/group_all_Ks.mat']);
group_parts=group_parts.affils;
% load in SNR masks
l_l = read_label([],'/cbica/projects/pinesParcels/data/H_SNR_masks/lh.Mask_SNR.label');
l_r = read_label([],'/cbica/projects/pinesParcels/data/H_SNR_masks/rh.Mask_SNR.label');
% +1 because matlab indexing starts differently than freesurfer
l_l_ind = l_l(:,1) + 1;
l_r_ind = l_r(:,1) + 1;
% check to make sure that mask indices match 0s in group consensus to ensure consistent masking throughout
% NOTE THAT THIS MASK FILE INDICATES THE PRESENCE OF VERTICES TO BE MASKED, NOT IN 0 = BAD 1 = GOOD FORMAT
if sum(group_parts(l_l_ind))~=0
	disp('you screwed up the left hemisphere mask numbnuts')
	exit(1);
else
end
if sum(group_parts(10242+l_r_ind))~=0
	disp('you screwed up the right hemisphere mask numbnuts')
	exit(1);
else
end
% change mask from 0s from 1 at shitty vertices to 0
surfMask.l = ones(10242,1);
surfMask.l(l_l_ind) = 0;
surfMask.r = ones(10242,1);
surfMask.r(l_r_ind) = 0;

% send off the group partitions for comparison to single-subject
group_parts=load([ProjectFolder '/SingleAtlas_Analysis/group_all_Ks.mat']);
group_parts=group_parts.affils;
group_parts_masked=group_parts(any(group_parts,2),:);

% only 97 new subjs
% test on 2
for s=1:97
	outdir = ['/cbica/projects/pinesParcels/data_psy/CombinedData/' num2str(subjs(s)) '/fc_metrics.mat']; 
	outdirp = ['/cbica/projects/pinesParcels/data_psy/CombinedData/' num2str(subjs(s)) '/pc_metrics.mat'];
	% commented out for overwriting
	%if ~exist(outdir, 'file')
		configfp=['/cbica/projects/pinesParcels/data_psy/CombinedData/' num2str(subjs(s)) '/fc_config.mat'];
		% save needed arguments
		save(configfp, 's', 'surfMask', 'Krange', 'subjs', 'group_parts_masked', 'outdir', 'outdirp');
		% turn this into a qsub command
		cmd = ['/cbica/projects/pinesParcels/multiscale_psy/scripts/run_subj_vert_fc_matlab_compile.sh $MATLAB_DIR ' configfp ];
		fid=fopen(['/cbica/projects/pinesParcels/data_psy/CombinedData/' num2str(subjs(s)) '/tmp.sh'], 'w');
		fprintf(fid,cmd);
		system(['qsub -l h_vmem=20G ' '/cbica/projects/pinesParcels/data_psy/CombinedData/' num2str(subjs(s)) '/tmp.sh']);
	% commented out for overwriting
	%end
end
