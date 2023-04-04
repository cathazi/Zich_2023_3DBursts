clear all;clc;close all;

% DEFINE PATHS
MAIN_PATH_1 = '/Volumes/My Book/UCL_MEG_P2/';
MAIN_PATH_2 = '/Users/catharinazich/Desktop/';

PATH_SCRIPT = [MAIN_PATH_2 'UCL_P2_45/scripts_TSS2/'];
PATH_DATA = [MAIN_PATH_1 'meg_data/'];
PATH_ROI = [MAIN_PATH_1 'AAL_MNI_ROIs/'];
PATH_FIG = [MAIN_PATH_2 'UCL_P2_45/figures/'];
PATH_STAT = [MAIN_PATH_2 'UCL_P2_45/stats/'];

SPM_PATH     = [MAIN_PATH_2 'Software/spm12/'];
OSL_PATH     = [MAIN_PATH_2 'Software/osl/'];
TOOL_PATH    = [MAIN_PATH_2 'UCL_P2_45/toolboxes/'];
addpath(genpath(PATH_SCRIPT));addpath(genpath(TOOL_PATH));
addpath(genpath(OSL_PATH));cd([OSL_PATH '/osl-core']);osl_startup;
rmpath(genpath([OSL_PATH 'spm12/external/fieldtrip/compat/']));
rmpath(genpath([OSL_PATH 'spm12/external/yokogawa_meg_reader/']));

% DEFINE JOBS
run.N_of_epochs  = 0;
run.roi          = 0;

run.preproc      = 0;
run.ica_id       = 0;

run.epoch        = 0;
run.beam         = 0;
run.tf           = 0;
run.cluster      = 0;
run.pca          = 0;
run.intTS        = 1;

% DEFINE PARAMS
SUBJ = {'8','3','7','1','6','2','4','5'};
time_epoch  = [-2000 2000]; % epoch boundaries
time_bl = [-1.8 -1.1]; % baseline for TF
freqs = [13:1:30];
thresh_std = [1:0.25:9];

cols_ss = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;153,153,153;166,86,40;247,129,191]./255;

%% get N of Epochs per Subj
if run.N_of_epochs
    for s = 1:length(SUBJ)
        N_of_epoch(s) = 0;
        PATH_SUBJ = [PATH_DATA SUBJ{s} '/New/'];cd(PATH_SUBJ);DAYS = dir('MEG_0*');
        for d = 1:length(DAYS)
            PATH_DAY = [PATH_SUBJ DAYS(d).name '/'];cd(PATH_DAY);SESS = dir('13-30_*.mat');
            for ss = 1:length(SESS)
                D = spm_eeg_load([PATH_DAY 'rtf_bre_' SESS(ss).name]);
                if length(D.badtrials) > 0
                    keyboard;
                end
                N_of_epoch(s) = N_of_epoch(s)+size(D,4);
            end %session
        end % day
    end % subj
    save([PATH_DATA 'N_of_epochs.mat'],'N_of_epoch');
end

%% 
for s = 1:length(SUBJ)
    PATH_SUBJ = [PATH_DATA SUBJ{s} '/New/'];mkdir(PATH_SUBJ);
    PATH_MRI = [PATH_SUBJ 'MRI/'];
    PATH_SURF = [PATH_SUBJ 'Surface/'];
    cd(PATH_SUBJ);DAYS = dir('MEG_0*');
    
    %% trasform Search Space ROI (SS ROI) and SensoriMotor ROI (SM ROI) and get SS ROI coords
    if run.roi
        % remove spm path and add spm path
        rmpath(genpath(OSL_PATH));
        addpath(genpath(SPM_PATH));
        % initialise spm and marsbar
        spm fmri
        marsbar
        
        cd(PATH_MRI);MRI = dir('*HeadCast*');PATH_MRI_SUBJ = [PATH_MRI MRI.name '/'];
        get_roi_coords(PATH_SCRIPT,PATH_MRI_SUBJ,PATH_ROI);
        % close spm, remove path and add osl path
        close all
        rmpath(genpath(SPM_PATH));
        addpath(genpath(OSL_PATH));cd([OSL_PATH '/osl-core']);osl_startup;
        rmpath(genpath([OSL_PATH '/spm12/external/fieldtrip/compat/']));
    end
    
    %% get and save all spatial infos
    [sl,sl_vi,SS_ROI_ind,ROI_ind] = get_spat_info(PATH_SUBJ,PATH_SURF,PATH_MRI,SUBJ{s});
    
    % reset defaults for cluster analysis
    if run.cluster
        cl_counter = 0;
        clear cl cl_info
        cl = cell(4000,1);
        cl(:) = {false(length(ROI_ind),18,72)};
    end
    
    for d = 1:length(DAYS)
        PATH_DAY = [PATH_SUBJ DAYS(d).name '/'];
        cd(PATH_DAY);
        if run.preproc
            SESS = dir('*.ds'); temp = [];for k=1:length(SESS); if SESS(k).isdir == 0 temp(k)=k; end; end; temp(temp==0)=[]; SESS(temp)=[];
        else
            SESS = dir('13-30_*.mat');
        end

        for ss = 1:length(SESS)
            if run.preproc
                %% IMPORT
                parts = strsplit(SESS(ss).name,{'_','.'});
                D = osl_import([PATH_DAY SESS(ss).name],'outfile',[PATH_DAY 'spm_' SUBJ{s} '_sess_0' num2str(d) '_run_' parts{4}]);
                
                %% COREGISTRATION
                cd(PATH_MRI_SUBJ);
                if strcmp(SUBJ{s},'gb')
                    T1 = dir('mm*.img');
                else
                    T1 = dir('s2*.nii');
                end
                
                if strcmp(SUBJ{s},'ad_no_helmet')
                    fids = get_fids;
                    matlabbatch{1}.spm.meeg.source.headmodel.D       = {D.fullfile};
                    matlabbatch{1}.spm.meeg.source.headmodel.val     = 1;
                    matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
                    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.mri = {[T1.folder '/' T1.name ',1']};
                    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres    = 3;
                    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
                    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
                    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
                    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = fids(s).nas;
                    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = fids(s).lpa;
                    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = fids(s).rpa;
                    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
                    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
                    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
                    spm_jobman('run', matlabbatch);
                end
                
                % check coregistration
                D = spm_eeg_load(D.fullfile);
                D.save;
                hf = check_coreg(D,fids(s));
                cd(PATH_MRI_SUBJ);
                savefig(hf,[SUBJ{s},'_sess_',num2str(d),'_block_',num2str(ss),'.fig']);
                close(hf);
                
                %% DOWNSAMPLE & FILTER BROAD BAND
                D = do_df(PATH_DAY,D.fname);
                
                %% ICA
                if D.montage('getnumber')>0
                    D = D.montage('remove',1:D.montage('getnumber'));
                end
                % detect bad segements
                D = osl_detect_artefacts(D,'badchannels',true,'badtime',true,'dummy_epoch_tsize',1);
                % run ica & automatic bad component detection
                D = osl_africa(D,'do_ica',true,'do_ident','auto','do_remove',true,'precompute_topos',true);
                
                % manual bad component detection
                if run.ica_id
                    D = D.montage('switch',0);
                    D = osl_africa(D,'do_ica',false','do_ident','manual','do_remove',true');
                end
                D.save
                
                %% FILTER BETA BAND
                D           = D.montage('switch',2);
                D.save
                S           = [];
                S.D         =  D;
                S.band      = 'bandpass';
                S.freq      = [13 30];
                S.dir       = 'twopass';
                S.prefix    = '13-30_';
                D = spm_eeg_filter(S);
            end % if preproc

            %% EPOCH
            if run.epoch
                epoch.sufix = 'e_';epoch.type = 'RT';epoch.value=66;
                D_cont = spm_eeg_load([PATH_DAY SESS(ss).name]);
                % add events
                [D_cont,eve_info,RT_temp] = add_trigger(D_cont,PATH_DAY);
                % epoch
                D = do_epoch(D_cont,time_epoch,epoch);
                % remove bad trials
                S_roi = [];S_roi.D = D;
                [D,goodtrials] = spm_eeg_remove_bad_trials(S_roi); % change org function so that goodtrials is a output
                % delete file previous file
                parts = strsplit(S_roi.D.fullfile,'.');delete(strcat(parts{1},'.*'));
                if size(D,3) ~= length(RT_temp(goodtrials))
                    keyboard;
                end
            end
            
            % these data are on OSF: https://doi.org/10.17605/OSF.IO/9CQHB

            %% BEAMFORM
            if run.beam
                D = spm_eeg_load([PATH_DAY 're_' SESS(ss).name]);
                % removes former beamformer montages
                if D.montage('getnumber')>1
                    D = D.montage('remove',2:D.montage('getnumber'));
                end
                D = D.montage('switch',1);
                D.save;
                
                S_roi = [];S_roi.D = D;
                D = osl_inverse_model(S_roi.D,sl.face_center(SS_ROI_ind,:),'prefix','b','pca_order',55);
                % delete file previous file
                parts = strsplit(S_roi.D.fullfile,'.');delete(strcat(parts{1},'.*'));
            end
            
            %% TF
            if run.tf
                D = spm_eeg_load([PATH_DAY 'bre_' SESS(ss).name]);
                D = D.montage('switch',D.montage('getnumber'));
                S_roi = []; S_roi.D = D;
                D = do_tf(D,freqs,time_bl);
            end
            
            %% 3D BUrsts CLUSTER
            % needs to be adapted if staptial domain is within ROI boarder or across (also within funciton select_clusters)
            if run.cluster
                vis_cl_check = 0;
                vis_cl_viewer = 0;
                
                D = spm_eeg_load([PATH_DAY 'rtf_bre_' SESS(ss).name]);
                tf_3d = squeeze(median(D(ROI_ind,:,:,:),1)); % avg tf over SM ROI (when indexing D no SS_ROI_ind necessary as D = SS_ROI_ind)
                parts = strsplit(SESS(ss).name,{'spm_','.mat'});
                
                % get threshold per session
                dat_temp = D(ROI_ind,:,:,:); % based on ROI only
                dat_pos = dat_temp > 0; % based on positive values only
                dat_temp = dat_temp(dat_pos);
                for t = 1:length(thresh_std)
                    thresh(t) = mean(dat_temp(:))+std(dat_temp(:))*thresh_std(t);
                    vpath_cS = D(ROI_ind,:,:,:) > thresh(t);
                    vpath_cS = squeeze(sum(vpath_cS,1));
                    [sim_corr{d}{ss}(t),h,p,ci] = pointbiserial(vpath_cS(:),tf_3d(:)); % correlation
                end
                thresh_idx = min(find(max(sim_corr{d}{ss})==sim_corr{d}{ss}));
                %thresh_idx = 9; % same threshold across sessions, days and subjects
                
                % get clusters per epoch
                for e = 1:size(D,4)
                    disp(['%%% Running cluster: ' parts{2} ', Epoch: ' num2str(e) ' %%%'])
                    tic
                    % threshold epoch
                    vpath_cS = D(ROI_ind,:,:,e) > thresh(thresh_idx);
                    %vpath_cS = D(:,:,:,e) > thresh(thresh_idx);
                    
                    % get f,t clusters
                    cl_ft = bwlabel(squeeze(sum(vpath_cS,1)));
                    %figure;subplot(121);imagesc(cl_ft);colorbar;subplot(122);imagesc(squeeze(mean(D(:,:,:,e),1)));colorbar;
                    
                    for c = 1:max(cl_ft,[],'all') % loop through t,f clusters
                        % get indces of current t,f cluster
                        [f_ind,t_ind] = ind2sub(size(cl_ft),find(c == cl_ft));
                        
                        % get surface cl for each cell that is part of current t,f cl
                        cl_surf = zeros(size(vpath_cS));
                        for z = 1:length(f_ind)
                            % cl_surf = zeros(length(SS_ROI_ind),1);
                            ptCloud = pointCloud(sl.face_center(SS_ROI_ind(ROI_ind(vpath_cS(:,f_ind(z),t_ind(z)))),:)); % limited to ROI
                            %ptCloud = pointCloud([sl.face_center(SS_ROI_ind(vpath_cS(:,f_ind(z),t_ind(z))),:)]); % outside ROI
                            [~, ptDist] = knnsearch(sl.face_center(SS_ROI_ind(ROI_ind)),sl.face_center(SS_ROI_ind(ROI_ind)),'K',2,'Distance','euclidean');
                            [labels,~] = pcsegdist(ptCloud,8);%max(ptDist(:,2));2:10
                            cl_surf(vpath_cS(:,f_ind(z),t_ind(z)),f_ind(z),t_ind(z)) = labels+f_ind(z)*10000+t_ind(z)*100;
                        end
                        if vis_cl_check
                            cl_surf_org = cl_surf;
                        end
                        
                        % combine surf clusters
                        cl_surf = combine_surf_clusters(sl.face_center,SS_ROI_ind(ROI_ind),cl_surf,t_ind,f_ind); % limited to ROI
                        %cl_surf = combine_surf_clusters(sl.face_center,SS_ROI_ind,cl_surf,t_ind,f_ind); % outside ROI
                        
                        % save each cl in one cell
                        [N,cl_indi_ind] = histcounts(cl_surf,[0:1:max(cl_surf,[],'all')]);
                        cl_indi_ind(1) = [];N(1)=[]; % remove zeros - not clusters
                        cl_indi_ind(N<2^3)=[]; % remove clusters that are too small part 1 (min 2 in each dim == min 8)
                        
                        for i = 1:length(cl_indi_ind)
                            % check if cl is big enough in each dimension and does not
                            % exceede the borders in any dimension
                            cl_indi = select_clusters(sl.face_center,SS_ROI_ind,ROI_ind,vpath_cS,cl_surf==cl_indi_ind(i),2,cl_ft,D,e); % change depending on limited to or outside ROI!
                            if ~isnan(cl_indi(1,1,1))
                                cl_counter = cl_counter+1;
                                cl{cl_counter} = cl_indi>0; % much faster when in cell array compared to 4D
                                % get mains
                                cl_info.epoch_index(cl_counter) = e;
                                cl_info.unique_index(cl_counter) = str2num([num2str(d) num2str(ss) sprintf('%02d',e)]);
                                cl_info.thresh(cl_counter) = thresh(thresh_idx);
                                cl_info.thresh_std(cl_counter) = thresh_std(thresh_idx);
                                cl_info.sim(cl_counter) = max(sim_corr{d}{ss});
                                cl_info = get_cluster_mains(cl_info,cl_counter,cl_indi,D(ROI_ind,:,:,e),D,sl,SS_ROI_ind(ROI_ind)); % limited to ROI
                                %cl_info = get_cluster_mains(cl_info,cl_counter,cl_indi,D(:,:,:,e),D,sl,SS_ROI_ind); % outside ROI
                            end
                        end % for individual cl witin tf cluster
                    end % for cluster (c)
                    toc
                end % for epoch (e)
                
                % save
                cl(length(cl_info.a_mean)+1:end) = [];
                save([PATH_SUBJ 'clusters_M1.mat'],'cl_info','cl','thresh','sim_corr','-v7.3');
            end % if
        end % session
    end % day
end %Subj


%% PCAs
if run.pca
    % PCA 1: get 1st PC across Tdur,Fspread,Ssize,i.e. 'Extend'
    % PCA 2: get 1st PC acrosss x,y,z of face_center, i.e. S_center
    % get data
    for s = 1:length(SUBJ)
        PATH_SUBJ = [PATH_DATA SUBJ{s} '/New/'];
        load([PATH_SUBJ 'clusters_M1.mat'],'cl_info');
        load([PATH_SUBJ 'spat_info.mat']);

        % for pca 1 (Extend)
        for_pca_11{s}(1,:) = zscore(cl_info.t_dur);
        for_pca_11{s}(2,:) = zscore(cl_info.f_spread);
        for_pca_11{s}(3,:) = zscore(cl_info.s_size);
        %for_pca_11{s}(4,:) = zscore(cl_info.a_q95);

        % for pca 2 (Space)
        for c = 1:length(cl_info.a_mean)
            parts = strsplit(cl_info.source_file{c},'rtf_');
            if c == 1
                D = spm_eeg_load([parts{1} parts{2}]);
            else
                if ~strcmp(D.fullfile,[parts{1} parts{2}])
                    D = spm_eeg_load([parts{1} parts{2}]);
                end
            end
            % transform coord of face_center into MNI space
            temp = sl_vi.face_center(cl_info.s_center_face_ind(c),:);
            MNI_coords = D.inv{1}.mesh.Affine(1:3,1:3) * temp' + D.inv{1}.mesh.Affine(1:3,4);
            %MNI_coords = D.inv{1}.mesh.Affine(1:3,1:3) *
            %cl_info.s_center_coord(:,c) + D.inv{1}.mesh.Affine(1:3,4); %
            % non inflated surface
            for_pca_12{s}(:,c) = MNI_coords;
        end
    end
    
    % combine data across subjects & run pca 1
    for_pca_2 = horzcat(for_pca_11{:})';
    [coeff_1,score_1,~,~,explained_1,mu_1] = pca(zscore(for_pca_2(:,1:3),1,1));
    vbls = {'T','F','S'};figure;biplot(double(coeff_1(:,1:3)),'VarLabels',vbls);
    for_pca_31 = score_1(:,1);
    
    % combine data across subjects & run pca 2
    for_pca_2 = horzcat(for_pca_12{:})';
    [coeff_2,score_2,~,~,explained_2,mu_2] = pca(zscore(for_pca_2,1,1));
    vbls = {'X','Y','Z'};biplot(double(coeff_2(:,1:3)),'VarLabels',vbls);
    for_pca_32 = score_2(:,1);
    for_pca_33 = score_2(:,2);
    for_pca_34 = score_2(:,3);
    save([PATH_DATA 'pcas'],'coeff_1','score_1','explained_1','mu_1','coeff_2','score_2','explained_2','mu_2');
    
    % split in subjects/clusters
    counter = 0;
    for s = 1:length(SUBJ)
        PATH_SUBJ = [PATH_DATA SUBJ{s} '/New/'];
        load([PATH_SUBJ 'clusters_M1.mat']);
        %cl_info = rmfield(cl_info,'s_center_face_pc');
        for_pca_41 = for_pca_31(counter+1:counter+length(for_pca_11{s}));
        for_pca_42 = for_pca_32(counter+1:counter+length(for_pca_11{s}));
        for_pca_43 = for_pca_33(counter+1:counter+length(for_pca_11{s}));
        for_pca_44 = for_pca_34(counter+1:counter+length(for_pca_11{s}));
        counter = counter + length(for_pca_11{s});
        for c= 1:length(for_pca_41)
            cl_info.extend(c) = for_pca_41(c);
            cl_info.s_center_face_pc1(c) = for_pca_42(c);
            cl_info.s_center_face_pc2(c) = for_pca_43(c);
            cl_info.s_center_face_pc3(c) = for_pca_44(c);
        end
        save([PATH_SUBJ 'clusters_M1.mat'],'cl_info','cl','thresh','sim_corr','-v7.3')
    end
end


