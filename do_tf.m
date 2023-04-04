function [D,tf_blc] = do_tf(D,freqs,time_bl)

% normalise (in MEG-TMS done in function where bad emg trials removed)
D = D.montage('switch',D.montage('getnumber'));
if D.montage('getnumber')>0
    % apply z-transformation (thus basleine correction is not only done for the mean )
    dat=D(:,:,:);
    dat=zscore(dat,[],[2 3]);
    outfile = fullfile(D.path,D.fname);
    Dnode = clone(montage(D,'switch',0),outfile,[size(dat,1),D.nsamples,size(dat,3)]);
    Dnode = chantype(Dnode,1:Dnode.nchannels,'VE');
    Dnode(:,:,:)=dat;
    D = Dnode;
    D.save
end

% TF
S = [];
S.D = D;
S.frequencies = freqs;
S.timewin = [-Inf Inf];
S.phase = 0;
S.method = 'mtmconvol';
S.settings.taper = 'dpss';
S.settings.timeres = 400;
S.settings.timestep = 50;
S.settings.freqres = 2.5;
D = spm_eeg_tf(S);
% delete previous files
% parts = strsplit(S.D.fullfile,'.');delete(strcat(parts{1},'.*'));
% NEED THIS FILE FOR PHASE ANALYSIS

% blc
S = [];
S.D = D;
S.method = 'Diff';
S.timewin = [time_bl(1)*1000 time_bl(2)*1000];
S.pooledbaseline = 0;
D = spm_eeg_tf_rescale(S);
% delete previous files 
parts = strsplit(S.D.fullfile,'.');delete(strcat(parts{1},'.*'));