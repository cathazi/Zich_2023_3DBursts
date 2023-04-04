function [TH_var,TH_mean,Rsq_var,Rsq_mean,Time,Speed_sl,Speed_sl_vi,Distance_sl] = intTS(cl,cl_info,sl_vi,sl,SS_ROI_ind,ROI_ind,sl_dist,sl_vi_dist)


temp = 1:size(cl,1);

for c = 1:length(temp)
    tic
    % for trouble shooting
    disp(['%%% Running cl: ',num2str(c)]);
    c = temp(c);

    parts = strsplit(cl_info.source_file{c},'rtf_');
    if c == temp(1) || ~strcmp(parts{2},D.fname)
        % load correct data set
        D = spm_eeg_load([parts{1} parts{2}]);

        % get envelope threshold (one per dataset)
        D_temp = squeeze(mean(D(ROI_ind,:,:),1));D_temp = reshape(D_temp,1,[]);D_temp_eve = envelope(D_temp);
        thresh_eve = mean(D_temp_eve);
    end

    % where in space is the burst active? (across all times and space) FORMAT: index
    spatial_mask = sum(cl{c},[2 3])>0;
    spatial_mask = spatial_mask.*[1:length(spatial_mask)]';spatial_mask(spatial_mask==0)=[];

    % when in time is the burst active? (specifically get first and last time index)
    [~,time_on] = min(abs(D.time-cl_info.t_onset(c)));
    [~,time_off] = min(abs(D.time-cl_info.t_offset(c)));

    if cl_info.epoch_index(c)
        % get data (prior to tf analysis) for that burst
        data = D(ROI_ind(spatial_mask),time_on:time_off,cl_info.epoch_index(c)); % limites to ROI
        %data = D(spatial_mask,time_on:time_off,cl_info.epoch_index(c)); % outside ROI

        % check and correct sign ambiguity
        data = check_sign_ambiguity(data,c,0);

        [TH_var_temp,TH_mean_temp,Rsq_var_temp,Rsq_mean_temp,Time_temp,Speed_sl_temp,Speed_sl_vi_temp,Distance_sl_temp] = get_tw_dir(data,sl_vi,sl,sl_dist,sl_vi_dist,SS_ROI_ind,ROI_ind,spatial_mask,D.inv{1}.mesh.Affine,thresh_eve);
        TH_var{c} = TH_var_temp;
        TH_mean{c} = TH_mean_temp;
        Rsq_var{c} = Rsq_var_temp;
        Rsq_mean{c} = Rsq_mean_temp;
        Time{c} = Time_temp;
        Speed_sl{c} = Speed_sl_temp;
        Speed_sl_vi{c} = Speed_sl_vi_temp;
        Distance_sl{c} = Distance_sl_temp;
    end
    toc
end





