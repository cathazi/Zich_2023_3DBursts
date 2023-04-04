function cl_info = get_cluster_mains(cl_info,cl_counter,cl_indi,envelope_cS,D,sl,SS_ROI_ind_temp)

cl_info.source_file{cl_counter} = D.fullfile;

% AMPLITUDE
cl_info.a_mean(cl_counter) = mean(envelope_cS(cl_indi>0),'all');
cl_info.a_q95(cl_counter) = quantile(reshape(envelope_cS(cl_indi>0),1,[]),0.95);

% TEMPORAL
temp = squeeze(sum(cl_indi,[1 2])>0).*[1:size(cl_indi,3)]'; temp(temp==0)=[];
cl_info.t_onset(cl_counter) = D.time(min(temp));
cl_info.t_offset(cl_counter) = D.time(max(temp));
cl_info.t_dur(cl_counter) = length(temp)*(D.time(2)-D.time(1));
cl_info.t_center(cl_counter) = cl_info.t_onset(cl_counter)+cl_info.t_dur(cl_counter)/2;
% interval time can not be computed for each burst individually

% SPECTRAL
temp = squeeze(sum(cl_indi,[1 3])>0).*[1:size(cl_indi,2)]; temp(temp==0)=[];
cl_info.f_Blow(cl_counter) = D.frequencies(min(temp));
cl_info.f_Bhigh(cl_counter) = D.frequencies(max(temp));
cl_info.f_spread(cl_counter) = length(temp)*(D.frequencies(2)-D.frequencies(1));
cl_info.f_center(cl_counter) = cl_info.f_Blow(cl_counter)+cl_info.f_spread(cl_counter)/2;

% SPATIAL
% total surface area
temp = zeros(1,size(sl.face_center,1));temp(SS_ROI_ind_temp) = squeeze(sum(cl_indi,[2 3])>0);temp = logical(temp);
cl_info.s_size(cl_counter) = get_areaIsosurface(sl.faces(temp,:),sl.vertices)*0.01; %unit: cm3

% x,y,z dimension (bounding box)
    %) get min-max in each dimension
x_ind = [find(min(sl.face_center(temp,1))==sl.face_center(:,1),1,'first') find(max(sl.face_center(temp,1))==sl.face_center(:,1),1,'first')];
y_ind = [find(min(sl.face_center(temp,2))==sl.face_center(:,2),1,'first') find(max(sl.face_center(temp,2))==sl.face_center(:,2),1,'first')];
z_ind = [find(min(sl.face_center(temp,3))==sl.face_center(:,3),1,'first') find(max(sl.face_center(temp,3))==sl.face_center(:,3),1,'first')];
    %) difference between bounding box coords
cl_info.s_size_x(cl_counter) = (max(sl.face_center([x_ind y_ind z_ind],1))-min(sl.face_center([x_ind y_ind z_ind],1)))*0.1; % unit cm
cl_info.s_size_y(cl_counter) = (max(sl.face_center([x_ind y_ind z_ind],2))-min(sl.face_center([x_ind y_ind z_ind],2)))*0.1;
cl_info.s_size_z(cl_counter) = (max(sl.face_center([x_ind y_ind z_ind],3))-min(sl.face_center([x_ind y_ind z_ind],3)))*0.1;

% center
    % get centers of mass
C_wb = mean(sl.face_center,1); % center of mass whole brain 
C_cl = mean(sl.face_center(temp,:),1); % center of mass burst
    % get intersection of vector through both centers of mass and surface
inter_temp = TriangleRayIntersection(C_wb, C_cl-C_wb, sl.vertices(sl.faces(:,1),:), sl.vertices(sl.faces(:,2),:), sl.vertices(sl.faces(:,3),:));
    % if more than one intercept take the one that is closest to center of mass
inter_temp = inter_temp.*[1:1:length(inter_temp)]';inter_temp(inter_temp==0)=[];
if length(inter_temp)>1
    for i = 1:length(inter_temp)
        inter_dist(i) = pdist([mean(sl.vertices(sl.faces(inter_temp(i),:),:),1);C_cl],'euclidean');
    end
end


cl_info.s_center_face_ind(cl_counter) = inter_temp(min(find(min(inter_dist)==inter_dist)));
cl_info.s_center_coord(:,cl_counter) = mean(sl.vertices(sl.faces(cl_info.s_center_face_ind(cl_counter),:),:),1);


%% VISUALIZE SPATIAL MAINS
% temp = double(temp);temp(temp==0)=NaN;
% f1 = figure('units','normalized','outerposition',[0 0 0.65 0.65],'DefaultAxesFontSize',18);hold on;
% subplot(1,3,[1 2]);hold on;
% %h = scatter3(face_center(SS_ROI_ind,1),face_center(SS_ROI_ind,2),face_center(SS_ROI_ind,3),36,[0.8 0.8 0.8],'filled');
% %set(h, 'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha', 0.2);        
% 
% % plot burst
% p = patch('Faces',sl.faces,'Vertices',sl.vertices,'FaceVertexCData',temp','EdgeColor','none','FaceColor','flat');
% set(p,'Facealpha',0.5);
% 
% % plot center of mass and vector
% scatter3(C_wb(1),C_wb(2),C_wb(3),'ko','filled');
% scatter3(C_cl(1),C_cl(2),C_cl(3),'ro','filled');
% plot3([C_wb(1) C_cl(1)],[C_wb(2) C_cl(2)],[C_wb(3) C_cl(3)],'r','Linewidth',2);
% 
% % plot center face on surface
% inter_face = nan(size(sl.faces,1),1);inter_face(cl_info.s_center_face_ind(cl_counter)) = 2;
% p = patch('Faces',sl.faces,'Vertices',sl.vertices,'FaceVertexCData',inter_face,'EdgeColor','none','FaceColor','flat');
% colormap(f1,hot);caxis([1 3]);%colorbar;
% scatter3(cl_info.s_center_coord(1,cl_counter),cl_info.s_center_coord(2,cl_counter),cl_info.s_center_coord(3,cl_counter),'r*');
% 
% % plot bounding box
% rec_x = [min(face_center([x_ind y_ind z_ind],1)) max(face_center([x_ind y_ind z_ind],1))];
% rec_y = [min(face_center([x_ind y_ind z_ind],2)) max(face_center([x_ind y_ind z_ind],2))];
% rec_z = [min(face_center([x_ind y_ind z_ind],3)) max(face_center([x_ind y_ind z_ind],3))];
% 
% plot3(rec_x,[min(rec_y) min(rec_y)],[min(rec_z) min(rec_z)],'g');
% plot3([min(rec_x) min(rec_x)],rec_y,[min(rec_z) min(rec_z)],'g');
% plot3([min(rec_x) min(rec_x)],[min(rec_y) min(rec_y)],rec_z,'g');
% 
% plot3(rec_x,[max(rec_y) max(rec_y)],[max(rec_z) max(rec_z)],'g');
% plot3([max(rec_x) max(rec_x)],rec_y,[max(rec_z) max(rec_z)],'g');
% plot3([max(rec_x) max(rec_x)],[max(rec_y) max(rec_y)],rec_z,'g');
% 
% plot3(rec_x,[min(rec_y) min(rec_y)],[max(rec_z) max(rec_z)],'g');
% plot3([min(rec_x) min(rec_x)],rec_y,[max(rec_z) max(rec_z)],'g');
% plot3([min(rec_x) min(rec_x)],[max(rec_y) max(rec_y)],rec_z,'g');
% 
% plot3(rec_x,[max(rec_y) max(rec_y)],[min(rec_z) min(rec_z)],'g');
% plot3([max(rec_x) max(rec_x)],rec_y,[min(rec_z) min(rec_z)],'g');
% plot3([max(rec_x) max(rec_x)],[min(rec_y) min(rec_y)],rec_z,'g');
% 
% % decorate
% xlabel('x - 1');ylabel('y -2');zlabel('z -3');
% title(['Area = ',num2str(cl_info.s_size(cl_counter)) ' cm^2']);
% axis equal;

%keyboard;



