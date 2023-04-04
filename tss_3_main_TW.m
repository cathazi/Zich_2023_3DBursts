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

% DEFINE PARAMS
SUBJ = {'8','3','7','1','6','2','4','5'};

%% INTERACTIONS
% need to get distances first
% run get_distances (own script, not function) using geodesic toolbox
if run.intTS
    for s = 1:length(SUBJ)
        disp(['%%% Running Subj: ',SUBJ{s} ' %%%'])
        PATH_SUBJ = [PATH_DATA SUBJ{s} '/New/'];
        load([PATH_SUBJ 'spat_info.mat']);
        load([PATH_SUBJ 'clusters_M1.mat']);
        load([PATH_SUBJ 'dist.mat']);

        [TH_var,TH_mean,Rsq_var,Rsq_mean,Time,Speed_sl,Speed_sl_vi,Distance_sl] = intTS(cl,cl_info,sl_vi,sl,SS_ROI_ind,ROI_ind,sl_dist,sl_vi_dist);
        save([PATH_SUBJ 'th_phi'],'TH_var','TH_mean','Rsq_var','Rsq_mean','Time','Speed_sl','Speed_sl_vi','Distance_sl','-v7.3');
    end %subj
end