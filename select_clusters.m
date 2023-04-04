function cl_indi = select_clusters(face_center,SS_ROI_ind,ROI_ind,vpath_cS,cl_surf,min_size_per_dim,cl_ft,D,e)
%% [cl_indi] = select_clusters(sl,SS_ROI_ind,vpath_cS,cl_surf,min_size_per_dim)
%

% Function selecting clusters
%
% Inputs
%   face_center: Nx3


% Take only clusters with min size in all dimensions
[v,f,t] = ind2sub(size(cl_surf),find(cl_surf==1));
if length(unique(v)) >= min_size_per_dim && length(unique(f)) >= min_size_per_dim && length(unique(t)) >= min_size_per_dim

    % Take only clusters that are within the the time, frequency boundaries
%    if min(t)>1 && max(t)<size(vpath_cS,3) % && min(f)>1 && max(f)<size(vpath_cS,2)
        
        %% Spatial selection
        % check if cluster is in SM ROI 
%        keyboard;
        on = sum(cl_surf,[2 3])>0;
        in_all = ismember(face_center(SS_ROI_ind(ROI_ind),:),face_center(SS_ROI_ind(on),:),'rows');
        
        % check if center is in SM ROI (center of mass is good enough here)
        %in_cent = inhull(mean(face_center(SS_ROI_ind_temp(on),:),1),roi_coords.mni');
       
%         %VISUALIZATION
% %         real_in = zeros(size(on));
% %         counter =0;
% %         for o=1:length(on)
% %            if on(o)==1
% %                counter = counter+1;
% %                if in_all(counter)==1
% %                     real_in(o)=1;
% %                end
% %            end
% %         end
% %         real_in=real_in>0;
% 
%         %k = boundary(face_center(SS_ROI_ind_temp,1),face_center(SS_ROI_ind_temp,2),face_center(SS_ROI_ind_temp,3),1);
%         f1 = figure('units','normalized','outerposition',[0 0 1 1]);
%         subplot(2,4,[1 6]);hold on;
%         h1 = scatter3(face_center(SS_ROI_ind,1),face_center(SS_ROI_ind,2),face_center(SS_ROI_ind,3),50,'o','MarkerFacecolor',[0.8 0.8 0.8],'MarkerEdgecolor',[0.8 0.8 0.8]);
%         set(h1,'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha', 0.2);
%         h2 = scatter3(face_center(SS_ROI_ind(ROI_ind),1),face_center(SS_ROI_ind(ROI_ind),2),face_center(SS_ROI_ind(ROI_ind),3),50,'o','MarkerFacecolor',[0.8 0.8 0.8],'MarkerEdgecolor',[0.8 0.8 0.8]);
%         %trisurf(k,roi_coords.mni(1,:)',roi_coords.mni(2,:)',roi_coords.mni(3,:)','Edgecolor','none','Facecolor',[0.4 0.8 1],'FaceAlpha',0.2)
%         scatter3(face_center(SS_ROI_ind(ROI_ind(on)),1),face_center(SS_ROI_ind(ROI_ind(on)),2),face_center(SS_ROI_ind(ROI_ind(on)),3),100,'ro','filled');
%         %scatter3(face_center(SS_ROI_ind_temp(real_in),1),face_center(SS_ROI_ind_temp(real_in),2),face_center(SS_ROI_ind_temp(real_in),3),100,'bo','filled');
%         %if in_cent
%         %    scatter3(mean(face_center(SS_ROI_ind_temp(on),1)),mean(face_center(SS_ROI_ind_temp(on),2)),mean(face_center(SS_ROI_ind_temp(on),3)),200,'go','filled');   
%         %else
%         %    scatter3(mean(face_center(SS_ROI_ind_temp(on),1)),mean(face_center(SS_ROI_ind_temp(on),2)),mean(face_center(SS_ROI_ind_temp(on),3)),200,'ko','filled');
%         %end
%         axis equal
%         
%         subplot(2,4,3);
%         %on = on .* [1:length(on)]';on(on==0)=[];
%         imagesc(D.time,D.frequencies,squeeze(mean(D(ROI_ind(on),:,:,e),1)));
%         colorbar;axis tight;set(gca,'YDir','normal');
%         ylabel('Frequency (Hz)');xlabel('Time (s)');title('TF current burst');
%         
%         subplot(2,4,4);
%         on = zeros(length(D.frequencies),length(D.time));
%         on(squeeze(sum(cl_surf,1)>0)) = 1;
%         imagesc(D.time,D.frequencies,on);
%         colorbar;axis tight;set(gca,'YDir','normal');
%         ylabel('Frequency (Hz)');xlabel('Time (s)');title('TF current burst');
%         
%         
%         subplot(2,4,7);
%         imagesc(D.time,D.frequencies,squeeze(mean(D(ROI_ind,:,:,e),1)));
%         colorbar;axis tight;set(gca,'YDir','normal');
%         ylabel('Frequency (Hz)');xlabel('Time (s)');title('TF whole epoch');
%         
%         subplot(2,4,8);
%         imagesc(D.time,D.frequencies,cl_ft);colorbar;axis tight;set(gca,'YDir','normal');
%         ylabel('Frequency (Hz)');xlabel('Time (s)');title('TF whole epoch');
%         
%         set(findall(gcf,'-property','FontSize'),'FontSize',26)
%         
%         keyboard;
%         close (f1)
        
        % take only those that have at least one voxel in ROI (center of mass is in ROI)
        if sum(in_all)>0 % in_cent
            cl_indi = zeros(squeeze(size(vpath_cS)));
            cl_indi(1:size(cl_surf,1),1:size(cl_surf,2),1:size(cl_surf,3)) = cl_surf;
         end %if(ven with SM)
%    end % if (time,freq boundaries)
end % if(size)

if exist('cl_indi','var') == 0
    cl_indi = NaN;
end


    