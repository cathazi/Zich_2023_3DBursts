function D = ucl_45_preproc01_df(path,file)

% downsample 300 Hz
S               = [];
S.D             = [path file];
S.fsample_new   = 300;
S.prefix        = 'df_'; % set this now, subsequent filters just overwrite
D = spm_eeg_downsample(S);
%delete(S.D);

% filter band
S.D             = [path 'df_' file];
S.band          = 'bandpass';
S.freq          = [1 95];
S.dir           = 'twopass';
S.prefix        = ''; % Overwrite previous file
D = spm_eeg_filter( S );

%filter notch
S               = [];
S.D             = [path 'df_' file];
S.band          = 'stop';
S.freq          = [49 51];
S.dir           = 'twopass';
S.prefix        = '';
D = spm_eeg_filter( S );

% S               = [];
% S.D             = [path 'df_' file];
% S.band          = 'stop';
% S.freq          = [99 101];
% S.dir           = 'twopass';
% S.prefix        = '';
% D = spm_eeg_filter( S );

% S               = [];
% S.D             = [path 'df_' file];
% S.band          = 'stop';
% S.freq          = [149 151];
% S.dir           = 'twopass';
% S.prefix        = '';
% D = spm_eeg_filter( S );

% S               = [];
% S.D             = [path 'df_' file];
% S.band          = 'stop';
% S.freq          = [199 201];
% S.dir           = 'twopass';
% S.prefix        = '';
% D = spm_eeg_filter( S );

% % quality check
% for d=1:size(D(:,:,:),1)3
%     if strcmp(D.chantype{d},'MEGMAG')||strcmp(D.chantype{d},'MEGPLANAR')||strcmp(D.chantype{d},'MEGGRAD')
%         [pxx,f]=pwelch(D(d,:,:),D.fsample,D.fsample/2,D.fsample,D.fsample);
%         pxx_log(d,:)=squeeze(10*log10(pxx));        
%     end
% end
% 
% figure;hold on;
% plot(f,pxx_log)
% plot([50 50],[min(min(pxx_log)) max(max(pxx_log))]);
% grid on;
% xlabel('Frequency (Hz)');ylabel('Power');
% title([file],'interp','none');
% print('-dpng',[path 'quality_check/filter' file '.png']);




