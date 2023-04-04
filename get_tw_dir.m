
function [TH_var_temp,TH_mean_temp,Rsq_var_temp,Rsq_mean_temp,Time_temp,Speed_sl_temp,Speed_sl_vi_temp,Distance_sl_temp] = get_tw_dir(data,sl_vi,sl,sl_dist,sl_vi_dist,SS_ROI_ind,ROI_ind,spatial_mask,affine,thresh_eve)

%% prep
% interpolate
xq = [1:0.01:size(data,2)];
data_i = interp1([1:size(data,2)],data',xq,'spline');data_i = data_i';

% set mean to timeseries to zero
data_is = data_i - mean(data_i,'all');

%% Get control points of the average
data_av = mean(data_is,1);

% get peaks and troughs on avg
[~,p_locs] = findpeaks(smoothdata(data_av)); % need smoothdata for surrogate data, check if also ok in real data, or if too smooth
[~,t_locs] = findpeaks(smoothdata(data_av*-1));

% get apmplitude threshold for peaks & throughs
eve = envelope(data_av);
eve(1:100) = NaN;eve(end-100:end) = NaN;

% if conditions are met
if length(p_locs) > 0 && length(t_locs) > 0 && (sum(eve(p_locs)>thresh_eve)>1 || sum(eve(t_locs)>thresh_eve)>1)
    locs_av_info(:,1) = [p_locs,t_locs];% combine peaks and throughs
    locs_av_info(:,2) = [ones(1,length(p_locs)),ones(1,length(t_locs))*-1];% binary identification
    locs_av_info = sortrows(locs_av_info,1);
    
    % get mid points between peak and through and vice versa (= boundaries for single voxels)
    for i = 1:size(locs_av_info,1)-1
        locs_temp(i)=round(mean([locs_av_info(i,1),locs_av_info(i+1,1)]));
    end
    locs_temp(end+1) = 1;locs_temp(end+1) = size(data_is,2);locs_temp = sort(locs_temp);
    
    for i = 1:length(locs_av_info)
        locs_av_info(i,3) = locs_temp(i);
        locs_av_info(i,4) = locs_temp(i+1);
    end
    
    %     % control figure
    %     figure;hold on;
    %     plot(data_is');
    %     plot(data_av,'k','Linewidth',2);
    %     plot([1 length(data_av)],[0 0],'k');
    %     for c = [1 3 4]
    %         for p = 1:length(locs_av_info)
    %             plot([locs_av_info(p,c) locs_av_info(p,c)],[min(data_av) max(data_av)],'k');
    %         end
    %     end
    %     axis tight;
    
    %% Get control points for each voxel
    for v = 1:size(data,1)
        clear locs_ss_*
        % get peak and through
        for i = 1:size(locs_av_info,1)
            if locs_av_info(i,2) == 1
                [~,locs_temp] = max(data_is(v,locs_av_info(i,3):locs_av_info(i,4)));
            else
                [~,locs_temp] = min(data_is(v,locs_av_info(i,3):locs_av_info(i,4)));
            end
            locs_ss_1(i) = locs_temp+locs_av_info(i,3)-1;
        end
        
        % get mid points
        for i = 1:length(locs_ss_1)-1
            locs_ss_2(i) = round(mean([locs_ss_1(i),locs_ss_1(i+1)]));
        end
        locs_ss(v,:) = sort([locs_ss_1,locs_ss_2]);
    end %voxel
    
    %     % control figure
    %     for i = 1:size(locs_ss,2)
    %        figure('units','normalized','outerposition',[0 0 1 1]);
    %        subplot(1,4,[1 3]);hold on;
    %        plot(data_is');
    %        plot(data_av,'k','Linewidth',2);
    %
    %        for v = 1:size(data_is,1)
    %             plot([locs_ss(v,i) locs_ss(v,i)],[min(data_av) max(data_av)]);
    %        end
    %        axis tight;
    %
    %        subplot(144);hold on;
    %        histogram(locs_ss(:,i));
    %     end
    
    %%
    % normalise each to the equivalent point obtained from the average
    locs_av = unique(sort([locs_av_info(:,1)',locs_av_info(:,3)',locs_av_info(:,4)']));
    locs_av(1) = [];locs_av(end) = [];
    
    if size(locs_av,2) == size(locs_ss,2) % solve differently!
        locs_rel = locs_av-locs_ss;
        % threshold (i.e. only take points above envelope threshold)
        locs_rel = locs_rel(:,eve(locs_av)>thresh_eve);
        locs_av = locs_av(eve(locs_av)>thresh_eve); % for visualisation

        % get gradient
        for i = 1:size(locs_rel,2)
                    % separete linear regressions
                    for e = 1:3 %x,y,z
                        % get regression coefficient (i.e. slope)
                        [reg_coef,S] = polyfit(locs_rel(:,i),sl_vi.face_center(SS_ROI_ind(ROI_ind(spatial_mask)),e),1); % reg_coef1(t,e,1) = slope, reg_coef1(t,e,2) = intercept
                        pts(i,e,1) = 0; pts(i,e,2) = reg_coef(1);
                        % R-sqaure
                        Rsq(i,e) = 1 - (S.normr/norm(sl_vi.face_center(SS_ROI_ind(ROI_ind(spatial_mask)),e) - mean(sl_vi.face_center(SS_ROI_ind(ROI_ind(spatial_mask)),e))))^2;
                    end
                     % normalise points to MNI to allow grand average across
                    % subj.
                    pts_aff(i,:,:) = affine(1:3,1:3) * squeeze(pts(i,:,:)) + affine(1:3,4);
            
%             % multiple regression
%             [b,~,~,~,stats] = regress(locs_rel(:,i),[ones(size(locs_rel(:,i))),sl_vi.face_center(SS_ROI_ind(ROI_ind(spatial_mask)),1),sl_vi.face_center(SS_ROI_ind(ROI_ind(spatial_mask)),2)]);
%             pts(i,[1 2],1) = [0 0];
%             pts(i,[1 2],2) = [b(2) b(3)];
%             % normalise points to MNI to allow grand average across
%             % subj.
%             pts_aff(i,:,:) = affine(1:2,1:2) * squeeze(pts(i,:,:)) + affine(1:2,4);
%             Rsq(i,[1 2]) = [stats(1)];
            
            % get speed
            % get delta time
            [val_min,idx_min] = min(locs_rel(:,i));
            [val_max,idx_max] = max(locs_rel(:,i));
            delta_time = (abs(val_min-val_max)/100)/300;% samples to sec
            
            % get delta space
            % not inflated
            [a1,b1] = knnsearch(sl_dist.face_center(sl_dist.ROI_ind,:),sl.face_center(SS_ROI_ind(ROI_ind(spatial_mask([idx_min]))),:));
            [a2,b2] = knnsearch(sl_dist.face_center(sl_dist.ROI_ind,:),sl.face_center(SS_ROI_ind(ROI_ind(spatial_mask([idx_max]))),:));
            
            Distance_sl(i) = mean([sl_dist.dist_geo(a1,a2),sl_dist.dist_geo(a2,a1)])* 0.001; %mm to m
            Speed_sl(i) = Distance_sl(i)/delta_time;
            %speed_sl_euc(i) = (sl_dist.dist_euc(a1,a2)*0.001)/delta_time;
            
            % inflated
            [a1,b1] = knnsearch(sl_vi_dist.face_center(sl_vi_dist.ROI_ind,:),sl.face_center(SS_ROI_ind(ROI_ind(spatial_mask([idx_min]))),:));
            [a2,b2] = knnsearch(sl_vi_dist.face_center(sl_vi_dist.ROI_ind,:),sl.face_center(SS_ROI_ind(ROI_ind(spatial_mask([idx_max]))),:));
            
            temp = mean([sl_vi_dist.dist_geo(a1,a2),sl_vi_dist.dist_geo(a2,a1)])* 0.001; %mm to m
            Speed_sl_vi(i) = temp/delta_time; %unit m/s
            %speed.sl_vi_euc(i) = (sl_vi_dist.dist_euc(a1,a2)*0.001)/delta_time;
        end
        
        UVW = pts_aff(:,:,2)-pts_aff(:,:,1);
        Rsq = mean(Rsq(:,1:2),2);% azimuth TH, elevation PHI
        % focus only on TH in the following therefore average over x,y for Rsq only
        
        %% get number of directions
        % remove z direction
        UVW_t = UVW(:,1:2);
        % get length of 2D UVW
        clear UVW_n
        for t = 1:size(UVW,1)
            UVW_n(t) = norm(UVW_t(t,:));
        end
        % get 2D unit vetors
        UVW_u = (1./UVW_n)'.*UVW_t;
        
        % remove outliers
        UVW_in = isoutlier(UVW_u,'gesd');
        UVW_in = sum(~UVW_in,2)==2;
        UVW_u = UVW_u(UVW_in,:);
        UVW = UVW(UVW_in,:);
        Rsq = Rsq(UVW_in);
        locs_av = locs_av(UVW_in);
        locs_rel = locs_rel(:,UVW_in);
        Speed_sl = Speed_sl(UVW_in);
        Speed_sl_vi = Speed_sl_vi(UVW_in);
        Distance_sl = Distance_sl(UVW_in);
        
        %% cluster directions
        % Spectral Clustering on Similarity Matrix
        % Get distance between each pair of observations using the pdist and squareform functions with the default Euclidean distance metric
        dist_temp = pdist(UVW_u);
        dist = squareform(dist_temp);
        % Construct the similarity matrix from the pairwise distance and check that the similarity matrix is symmetric
        S = exp(-dist.^2);
        if ~issymmetric(S)
            TH_var_temp = NaN;
            TH_mean_temp = NaN;
            Rsq_var_temp = NaN;
            Rsq_mean_temp = NaN;
            Time_temp = NaN;
            Speed_sl_temp = NaN;
            Speed_sl_vi_temp = NaN;
            Distance_sl_temp = NaN;
        else
            % Limit the similarity values to 0.7
            S(S<0.7) = 0;
            % Create graph object from S
            G_eps = graph(S);
            
            % Identify the number of connected groups
            labels = conncomp(G_eps);
            labels_uni = unique(labels);
            
            % remove groups that are too small
            for y = 1:length(labels_uni)
                if sum(labels==labels_uni(y)) < length(labels)/(length(labels_uni)*2)
                    labels(labels==labels_uni(y)) = NaN;
                end
            end
            labels_uni = unique(labels);
            
            % sort based on when the direcetion occures in time
            for y = 1:length(labels_uni(~isnan(labels_uni)))
                y_time(y,1) = mean(locs_av(labels == labels_uni(y)));
                y_time(y,2) = y;
            end
            y_time_order = sortrows(y_time,1);
            
            % extract metrics per direction
            for y = 1:size(y_time_order,1)
                TH_temp = cart2sph(UVW_u(labels == y_time_order(y,2),1),UVW_u(labels == y_time_order(y,2),2),zeros(sum(labels == y_time_order(y,2)),1));
                TH_var_temp(y) = circ_std(TH_temp);
                TH_mean_temp(y) = circ_mean(TH_temp,UVW_n(labels == y_time_order(y,2))');
                Rsq_var_temp(y) = std(Rsq(labels == y_time_order(y,2)));
                Rsq_mean_temp(y) = mean(Rsq(labels == y_time_order(y,2)));
                Time_temp(y) = y_time_order(y,1);
                Speed_sl_temp(y) = mean(Speed_sl(labels == y_time_order(y,2)));
                Speed_sl_vi_temp(y) = mean(Speed_sl_vi(labels == y_time_order(y,2)));
                Distance_sl_temp(y) = mean(Distance_sl(labels == y_time_order(y,2)));
                %[TH_mean2{c}(y),PHI{c}(y)] = cart2sph(mean(UVW(labels == labels_uni(y),1)),mean(UVW(labels == labels_uni(y),2)),mean(UVW(labels == labels_uni(y),3)));
            end
            
            %     % control figure
            %     for i = 1:size(UVW,1)
            %         UVW_n(i) = norm(UVW(i,:));
            %     end
            %     [~,i_idx] = max(UVW_n);
            %     f1 = figure;
            %     subplot(2,2,1);hold on;
            %     plot(data_is');
            %
            %     subplot(2,2,3);hold on;
            %     daspect manual;pbaspect manual
            %     plot([0 0],[-max(abs(pts(:,:,2)),[],'all') max(abs(pts(:,:,2)),[],'all')],'color',[0.4 0.4 0.4]);
            %     plot([-max(abs(pts(:,:,2)),[],'all') max(abs(pts(:,:,2)),[],'all')],[0 0],'color',[0.4 0.4 0.4]);
            %     arrow3(zeros(size(pts,1),3),pts(:,:,2),'k0.5',0.5,0.5);
            %     arrow3([0 0 0],pts(i_idx,:,2),'c2',1,1);
            %     arrow3([0 0 0],mean(pts(:,:,2),1),'m2',1,1);
            %     axis tight;axis square;view(0,90);
            %     %xticks([]);yticks([]);zticks([]);
            %
            %     temp = nan(size(sl_vi.face_center,1),1);
            %     temp(SS_ROI_ind(ROI_ind(spatial_mask))) = locs_rel(:,i_idx);
            %     subplot(2,2,2);hold on;
            %     patch('Faces',sl_vi.faces,'Vertices',sl_vi.vertices,'FaceVertexCData',temp,'EdgeColor','none','FaceColor','flat');
            %     b = colorbar;b.Location = 'east';b.Label.String = 'Phase';
            %     view(0,90);
            %     util_resize_patch(sl_vi,SS_ROI_ind,ROI_ind(spatial_mask));
            %     xticks([]);yticks([]);zticks([]);
            %
            %     temp(SS_ROI_ind(ROI_ind(spatial_mask))) = mean(locs_rel,2);
            %     subplot(2,2,4);hold on;
            %     patch('Faces',sl_vi.faces,'Vertices',sl_vi.vertices,'FaceVertexCData',temp,'EdgeColor','none','FaceColor','flat');
            %     b = colorbar;b.Location = 'east';b.Label.String = 'Phase';
            %     view(0,90);
            %     util_resize_patch(sl_vi,SS_ROI_ind,ROI_ind(spatial_mask));
            %     xticks([]);yticks([]);zticks([]);
            %     close(f1);
        end
    else
        TH_var_temp = NaN;
        TH_mean_temp = NaN;
        Rsq_var_temp = NaN;
        Rsq_mean_temp = NaN;
        Time_temp = NaN;
        Speed_sl_temp = NaN;
        Speed_sl_vi_temp = NaN;
        Distance_sl_temp = NaN;
    end
else
    TH_var_temp = NaN;
    TH_mean_temp = NaN;
    Rsq_var_temp = NaN;
    Rsq_mean_temp = NaN;
    Time_temp = NaN;
    Speed_sl_temp = NaN;
    Speed_sl_vi_temp = NaN;
    Distance_sl_temp = NaN;
end