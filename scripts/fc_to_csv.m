outdir='/cbica/projects/pinesParcels/results_psy/aggregated_data/fc/'

% to become 2 to 30 when stuff finishes running someday
Krange=2:30

subjs=load('/cbica/projects/pinesParcels/data_psy/PsychSubjs.txt');

% load in the cell struct array frakenmatrices
ind_file=load('/cbica/projects/pinesParcels/results_psy/aggregated_data/ind_conmats_allscales_allsubjs.mat');
ind_feats=ind_file.ind_mats;

% calculate dataframe size
win_over_scales=zeros((length(Krange)*((min(Krange)+max(Krange))/2)),1);
bw_over_scales=[];
for K=Krange
	bw_over_scales=[bw_over_scales, zeros(1,(((K-1)*(K))/2))];
end

% get length of each, multiply by 3 for the 3 operationalizations of FC, and make an empty array of those combined +1 for subj ids
full_df_colnum=(length(bw_over_scales)+1+length(win_over_scales));
% length(subjs) +1 so theres a row for colnames
% looks like we can expect 16,357 columns for this dataframe if Krange=2:30

%% will have to change all of these to cells to accom. mixed strings and numeric
df=cell(length(subjs)+1,full_df_colnum);

rownames=strings(length(subjs)+1,1);
rownames(1)='subjects' ;
rownames(2:(length(subjs)+1))=num2str(subjs);
df(:,1)=cellstr(rownames);

% make a colnames row
colnames=strings(full_df_colnum+1);
colnames(1)='bblid';
% thirds names - legacy variable name - only using indiviualized fc "third"
thirdsnames=["ind"];
% initialize an empty df for each of three fc eval methods, +1 for colnames and subjnames
third_df=cell(length(subjs)+1,((full_df_colnum-1))+1,1);

% for each third
%for t=1:3;
% just individualized for now
t=1;

thirdname=thirdsnames(t);
pertinent_feat_name=strjoin([thirdname '_feats'],'');
feats=eval([pertinent_feat_name]);

%% starting with WITHIN %%
win_strings=strings(length(win_over_scales),1);
% make index of where features for each scale belong in string
Kind_w={};
% make a df just for within network connectivities: creating different DFs for each type of FC feature to make this whole script more modular in order to help with troubleshooting
df_win=cell(length(subjs)+1,length(win_over_scales));
% make an index of which places in feature vec align with which scale
for K=Krange
	K_start=((K-1)*(K))/2;
	K_end=(((K-1)*(K))/2)+K-1;
	Kind_w{K}=K_start:K_end;
end
for K=Krange
	Kind=Kind_w{K};
	for N=1:K
		curindex=Kind(N);
		win_strings(curindex)=strcat(thirdname, '_win_FC_scale', num2str(K), '_net', num2str(N));
	end
	df_win(1,Kind)=cellstr(win_strings(Kind));
	% extract matrices from this scale
	featurematrix=feats{K};
	for s=1:length(subjs);
		subjmat=featurematrix(:,:,s);
		subjwin=diag(subjmat);
		% s+1 because 1st row is colnames
		df_win(s+1,Kind)=num2cell(subjwin);
	end
end
	
%% BETWEEN Net FC %%
bw_strings=strings(length(bw_over_scales),1);
% between-FC dataframe for modular org of this df
df_bw=cell(length(subjs)+1,length(bw_over_scales));
% make index of which places in feature vec align with which scale
Kind_bw={};
% instead of st8 triangular numbers (almost str8, n-1 instead of n) for within, between needs summed num of triang nums from all prev. scales
K_start=1;
K_end=1;
Kind_bw{2}=K_start:K_end;
% start for 2nd scale, which is 3, just assign 1 to first scale
for K=(Krange(2):(length(Krange)+1))
	K_start=K_end+1;
	% -1 because range is inclusive of start and end
	K_end=K_start+((K)*(K-1)/2)-1;
	Kind_bw{K}=K_start:K_end;
end
% now that the index is made, populate df with colnames and values
for K=Krange
	Kind=Kind_bw{K};
	netlist=1:K;
	% set up a prev-net-index tracker so next net can fill in where previous left off
	prevnetind=0;
	% for all between network connectivities at this scale
	for N=netlist;
		% get networks from netlist that have not bee iterated over 
		unrecordednets=netlist(netlist>N);
		% This N's indices (FC with all unrecorded nets), +1 because indices are inclusive
		curnetind=((prevnetind+1):(prevnetind+length(unrecordednets)));
		% and where it corresponds to in greater K indices
		wheretofill=Kind(curnetind);
		% record b/w FC from curnet and unrecorded nets
		bwstringsN=[];
		for b=unrecordednets;
			bwstringsN=[bwstringsN, strcat(thirdname,'_bw_FC_scale',num2str(K),'_nets',num2str(N),'_and_',num2str(b))];
		end
	bw_strings(wheretofill)=bwstringsN;
	% update ending point of previous network indices
	prevnetind=max(curnetind);
	end
	df_bw(1,Kind)=cellstr(bw_strings(Kind));
	% on to actual values
	featurematrix=feats{K};
	% now add in subject values to the df using same indices but for subj row
	for s=1:length(subjs);
		subjmat=featurematrix(:,:,s);
		% extract vals for each network, parallel to above
		prevnetind=0;
		% reset to null
		subjbwvals=[];
		for N=netlist;
			unrecordednets=netlist(netlist>N);
			curnetind=((prevnetind+1):(prevnetind+length(unrecordednets)));
			% without wheretofill beccause we are not populating a global vector here
			bwval_forloop=[];
			for b=unrecordednets;
                	      	subjbwval=subjmat(N,b);
				bwval_forloop=[bwval_forloop,subjbwval];
			end	
			prevnetind=max(curnetind);
			subjbwvals(curnetind)=bwval_forloop;
		end
		df_bw(s+1,Kind)=num2cell(subjbwvals);
	end
end
% merge all fc features into dataframe of this third
%       |____________BW_CONS__________|
%COLNAMES|
% SUBJ_1|
% SUBJ_N|
% Merge all the fc metrics - save first column for subject names	
third_df(:,2:full_df_colnum,t)=[df_win df_bw];
% Make subjs col for all 1/3rd dfs
third_df(2:(length(subjs)+1),1,t)=num2cell(subjs);	
third_df(1,1,t)=cellstr('Subjects');

% combine FC feature-type dfs for this third
df=[third_df(:,:,1)];
% save in R friendly format
writetable(cell2table(df),strcat(outdir,'master_fcfeats_psy.csv'));

