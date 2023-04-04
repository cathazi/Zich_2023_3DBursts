function D = do_epoch (D_cont,time_epoch,epoch)

% define epoch parameter
D_cont    = D_cont.montage('switch',D_cont.montage('getnumber'));
S.D       = D_cont;
S.timewin = [time_epoch(1) time_epoch(2)];

S.trialdef(1).conditionlabel = epoch.type;
S.trialdef(1).eventtype = epoch.type;
S.trialdef(1).eventvalue = epoch.value;
S.reviewtrials = 0;
S.save = 0;
S.epochinfo.padding = 0;
S.event = D_cont.events;
S.fsample = D_cont.fsample;
S.timeonset = D_cont.timeonset;
[epochinfo.trl,epochinfo.conditionlabels,~] = spm_eeg_definetrial(S);

% do epoching
S2        = epochinfo;
S2.D      = D_cont;
S2.prefix = epoch.sufix;
D         = osl_epoch(S2);