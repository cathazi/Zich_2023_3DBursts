% get distances between vertices (start seed all center faces in ROI)

clear all;clc;

global geodesic_library;
geodesic_library = 'geodesic_debug';      %"release" is faster and "debug" does additional checks
addpath(genpath('/Users/catharinazich/Desktop/UCL_P2_45/toolboxes/geodesic_matlab-master'));

MAIN_PATH_1 = '/Volumes/My Book/UCL_MEG_P2/';
PATH_DATA = [MAIN_PATH_1 'meg_data/'];
SUBJ = {'8','3','7','1','6','2','4','5'};

%% find x coordinate that separates the hemispheres
%only do once

hemi_sep_x = [8;-4;-7;19.5;0;-5;-5;-3];
% for s = 1:length(SUBJ)
%     PATH_SUBJ = [PATH_DATA SUBJ{s} '/New/'];
%     load([PATH_SUBJ 'spat_info.mat']);
% 
%     figure;hold on;
%     scatter3(sl_vi.vertices(:,1),sl_vi.vertices(:,2),sl_vi.vertices(:,3),'k*');
%     axis equal;view(0,90);
%     plot3([hemi_sep_x(s) hemi_sep_x(s)],[-100 120],[100 100],'g')
% end

%%
for s = 1:length(SUBJ)
    PATH_SUBJ = [PATH_DATA SUBJ{s} '/New/'];
    load([PATH_SUBJ 'spat_info.mat']);
    %% remove one hemisphere (as only one mesh at a time and the two hemis are 2 meshes)
    % (hemisphere speration needs to be done on inflated surf, as not easy
    % detectable on non-inflated surface)
    
    % make new indicies for verticies for left hemisphere only
    clear left_hemi_v_idx
    left_hemi_v_idx(1,:) = single(sl_vi.vertices(:,1) < hemi_sep_x(s));
    left_hemi_v_idx(2,:) = left_hemi_v_idx .* [1:length(left_hemi_v_idx)];
    counter = 0;
    for v = 1:size(left_hemi_v_idx,2)
        if left_hemi_v_idx(1,v) == 1
            counter = counter +1;
            left_hemi_v_idx(3,v) = counter;
        end
    end
    
    % change index of verticies in faces (these indices are the same for inflated and uninflated surface)
    clear sl_dist
    for f = 1:size(sl_vi.faces,1)
        for d = 1:3
            [~,b] = ismember(double(sl_vi.faces(f,d)),left_hemi_v_idx(2,:));
            if b > 0
                sl_dist.faces(f,d) = left_hemi_v_idx(3,b);
            end
        end
    end
    sl_dist.faces = sl_dist.faces(any(sl_dist.faces,2),:);
    
    % loop through surfaces (non inflated, inflated)
    for g = 1:2
        % get verticies
        if g == 1
            % for inflated surface
            sl_dist.vertices = sl_vi.vertices(logical(left_hemi_v_idx(1,:)),:);
            sl_xx = sl_vi;
        else
            % for non inflated surface
            sl_dist = rmfield(sl_dist,{'vertices','face_center','dist_geo','dist_euc','ROI_ind'});
            sl_dist.vertices = sl.vertices(logical(left_hemi_v_idx(1,:)),:);
            sl_xx = sl;
        end
        
        % get center of faces
        for f = 1:size(sl_dist.faces,1)
            sl_dist.face_center(f,:) = mean(sl_dist.vertices(sl_dist.faces(f,:),:));
        end
        
        % new ROI index as SS_ROI_ind comprsies more than only the left hemisphere
        % voxels
        sl_dist.ROI_ind = ismember(sl_dist.face_center,sl_xx.face_center(SS_ROI_ind(ROI_ind),:),'row');
        sl_dist.ROI_ind = sl_dist.ROI_ind.*[1:length(sl_dist.ROI_ind)]';sl_dist.ROI_ind(sl_dist.ROI_ind==0)=[];
        
        % control figure
        figure;hold on;
        scatter3(sl_dist.vertices(:,1),sl_dist.vertices(:,2),sl_dist.vertices(:,3),'k*');
        scatter3(sl_dist.face_center(sl_dist.ROI_ind,1),sl_dist.face_center(sl_dist.ROI_ind,2),sl_dist.face_center(sl_dist.ROI_ind,3),'r*');
        axis equal;view(0,90);
        
        % initialise mesh and algorithm
        mesh = geodesic_new_mesh(sl_dist.vertices,sl_dist.faces);
        algorithm = geodesic_new_algorithm(mesh, 'exact');
        
        %% get distance for each face_center in ROI to each face_center
        % for each face_center
        clear dist*
        for face_id1 = 1:length(sl_dist.ROI_ind)
            disp(['%%% Running ',num2str(face_id1),'/',num2str(length(sl_dist.ROI_ind)),'%%%'])
            
            % seed
            coords = sl_dist.face_center(sl_dist.ROI_ind(face_id1),:);
            source_points = {geodesic_create_surface_point('face',sl_dist.ROI_ind(face_id1),coords)};
            
            % propagation stage of the algorithm (the most time-consuming)
            geodesic_propagate(algorithm, source_points);
            
            % find distances to all vertices of the mesh; in this example we have a single source, so source_id is always equal to 1
            [~, distances] = geodesic_distance_and_source(algorithm);
            
            % control figure
            % figure;hold on;
            % p = patch('Faces',sl_dist.faces,'Vertices',sl_dist.vertices,'FaceVertexCData',distances,'EdgeColor','k','FaceColor','flat');
            % axis equal;view(0,90);colorbar;colormap('jet');caxis([0 70]);
            % c = scatter3(coords(1),coords(2),coords(3),100,'mo','filled');
            
            % distance to other faces in ROI (mean of verticies)
            parfor face_id2 = 1:length(sl_dist.ROI_ind)
                vertex_idx_temp = sl_dist.faces(sl_dist.ROI_ind(face_id2),:);
                
                if face_id1 == face_id2
                    dist_geo(face_id1,face_id2) = 0; % otherwise mean of distance between center of face to each vertex of face ~= 0
                else
                    dist_geo(face_id1,face_id2) = mean(distances(vertex_idx_temp));
                end
                dist_euc(face_id1,face_id2) = pdist([coords;mean(sl_dist.vertices(vertex_idx_temp,:),1)]);
                
                %v = scatter3(sl_dist.vertices(verts,1),sl_dist.vertices(verts,2),sl_dist.vertices(verts,3),50,'mo','filled');
                %destination = geodesic_create_surface_point('vertex',verts(1),sl_dist.vertices(verts(1),:));
                %path = geodesic_trace_back(algorithm, destination);     %find a shortest path from source to destination
                %[x,y,z] = extract_coordinates_from_path(path);                                  %prepare path data for plotting
                %h = plot3(x*1.001,y*1.001,z*1.001,'w-','LineWidth',2);    %plot path
            end % face_id2
        end % face_id1
        
        sl_dist.dist_geo = dist_geo; sl_dist.dist_euc = dist_euc; % hand over to srtucture outside parfor, rather than in parfor
        if g==1
            sl_vi_dist = sl_dist;
        else
            sl_dist = sl_dist;
        end
        
        geodesic_delete;
    end % inflated/ not inflated
    
    save([PATH_SUBJ 'dist'],'sl_dist','sl_vi_dist','-v7.3');
end  % subj




