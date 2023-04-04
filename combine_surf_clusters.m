function cl_surf = combine_surf_clusters(face_center,SS_ROI_ind,cl_surf,t_ind_org,f_ind_org)
vis_cl = 0;

%f_ind = flip(f_ind);t_ind = flip(t_ind);
for d = 1:4
   
    % check directions again!
    if d == 1
        [t_ind,idx] = sort(t_ind_org);
        f_ind = f_ind_org(idx);
    elseif d==2
        [t_ind,idx] = sort(t_ind_org);t_ind = flip(t_ind);
        f_ind = f_ind_org(idx);f_ind = flip(f_ind);
    elseif d==3
        [f_ind,idx] = sort(f_ind_org);
        t_ind = t_ind_org(idx);
    elseif d==4
        [f_ind,idx] = sort(f_ind_org);f_ind = flip(f_ind);
        t_ind = t_ind_org(idx);t_ind = flip(t_ind);
    end
    
    % combine surf clusters
    for z = 1:length(f_ind)
        % surface clusters in this t,f cell
        curr_data = cl_surf(:,f_ind(z),t_ind(z));
        % number of surface clusters in this t,f cell
        curr_cl_tot = unique(curr_data);curr_cl_tot(curr_cl_tot==0)=[];
        
        % vis
        if vis_cl
            cols = [31,120,180;51,160,44;227,26,28;255,127,0;202,178,214;106,61,154;255,255,153;177,89,40;166,206,227;178,223,138;251,154,153;253,191,111;...
            31,120,180;51,160,44;227,26,28;255,127,0;202,178,214;106,61,154;255,255,153;177,89,40;166,206,227;178,223,138;251,154,153;253,191,111]./255;

            f1 = figure('units','normalized','outerposition',[0 0 1 1],'DefaultAxesFontSize',18);
            subplot(3,3,5);hold on;
            for u=1:length(curr_cl_tot)
                scatter3(face_center(SS_ROI_ind(curr_data==curr_cl_tot(u)),1),face_center(SS_ROI_ind(curr_data==curr_cl_tot(u)),2),face_center(SS_ROI_ind(curr_data==curr_cl_tot(u)),3),36,cols(u,:),'filled');
            end
            grid on;axis square;view(-39,27);
            xlim([-80 80]);ylim([-100 100]);zlim([-60 100]);legend(num2str(curr_cl_tot));
            title(['t:' num2str(t_ind(z)) ', f:' num2str(f_ind(z))]);
        end
        
        conn = 8;
        [neighs,neighs_sb] = get_neighs(conn,t_ind,f_ind,z);
        
        for n = 1:size(neighs,1)
            % check if neighbours in bounds of whole t,f matrix
            % & if neigh is in current t,f cluster
            if (neighs(n,1)>=1 && neighs(n,1)<=size(cl_surf,2) && neighs(n,2)>=1 && neighs(n,2)<=size(cl_surf,3)) && ismember(neighs(n,:),[f_ind,t_ind],'rows')
                % surface clusters in neighbour t,f cell
                neigh_data = cl_surf(:,neighs(n,1),neighs(n,2));
                % number of surface clusters in neighbour t,f cell
                neigh_cl_tot = unique(neigh_data);neigh_cl_tot(neigh_cl_tot==0)=[];
                
                % vis
                if vis_cl
                    subplot(3,3,neighs_sb(n));hold on;
                    for u=1:length(neigh_cl_tot)
                        scatter3(face_center(SS_ROI_ind(neigh_data==neigh_cl_tot(u)),1),face_center(SS_ROI_ind(neigh_data==neigh_cl_tot(u)),2),face_center(SS_ROI_ind(neigh_data==neigh_cl_tot(u)),3),36,cols(u,:),'filled');
                    end
                    grid on;axis square;view(-39,27);
                    xlim([-80 80]);ylim([-100 100]);zlim([-60 100]);zlim([0 90]);legend(num2str(neigh_cl_tot));
                    title(['t:' num2str(neighs(n,2)) ', f:' num2str(neighs(n,1))]);
                end
                
                for curr_cl = 1:length(curr_cl_tot)
                    curr_mask = curr_data == curr_cl_tot(curr_cl);
                    for neigh_cl = 1:length(neigh_cl_tot)
                        neigh_mask = neigh_data == neigh_cl_tot(neigh_cl);
                        
                        if sum((curr_mask+neigh_mask)>1)>1 % if current & next mask have at least one connectioning point
                            
                            % vis
                            if vis_cl
                                f2 = figure('units','normalized','outerposition',[0 0 0.4 0.65],'DefaultAxesFontSize',18);hold on;
                                scatter3(face_center(SS_ROI_ind(curr_mask),1),face_center(SS_ROI_ind(curr_mask),2),face_center(SS_ROI_ind(curr_mask),3),36,cols(curr_cl,:),'filled');
                                scatter3(face_center(SS_ROI_ind(neigh_mask),1),face_center(SS_ROI_ind(neigh_mask),2),face_center(SS_ROI_ind(neigh_mask),3),36,cols(neigh_cl,:),'filled');
                                scatter3(face_center(SS_ROI_ind((curr_mask+neigh_mask)>1),1),face_center(SS_ROI_ind((curr_mask+neigh_mask)>1),2),face_center(SS_ROI_ind((curr_mask+neigh_mask)>1),3),36,'k','filled');
                                grid on;axis square;view(-39,27);
                                xlim([-80 80]);ylim([-100 100]);zlim([-60 100]);legend(['curr ',num2str(curr_cl_tot(curr_cl))],['neigh ',num2str(neigh_cl_tot(neigh_cl)) ],'ven');
                                close(f2)
                            end
                            
                            % relabel
                            neigh_data(neigh_mask) = curr_cl_tot(curr_cl);
                            cl_surf(neigh_mask,neighs(n,1),neighs(n,2)) = curr_cl_tot(curr_cl);
                        end
                    end % for cl neighbouring t,f cell (neigh_cl)
                end % for cl current t,f cell (curr_cl)
            end % if (neigh validity check)
        end % for neigh (n)
        
        if vis_cl
            close(f1)
        end
    end % for cell (z)
    
end % for directions (d)
