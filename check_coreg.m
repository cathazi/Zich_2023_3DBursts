function hf = check_coreg(D,fids)

hf = figure;
if any(strcmp(get(get(hf,'children'),'type'),'axes'))
    [az,el] = view;
else
    az = 125; el = 15;
end
ha = axes('parent',hf); hold on
rotate3d(hf,'ON');

cols = {'r','g','b'};
D = montage(D,'switch');
mesh = spm_eeg_inv_transform_mesh(D.inv{1}.datareg(1).fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh);
%mesh_xxx = spm_eeg_inv_transform_mesh(D.inv{1}.datareg(1).fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh.tess_mni);
%mesh_xxx = struct('faces',mesh_xxx.face,'vertices',mesh_xxx.vert);

fid_polhemus = D.inv{1}.datareg(1).fid_eeg.fid.pnt;
fid_mri = D.inv{1}.datareg(1).fid_mri.fid.pnt;
%patch(mesh_xxx,'FaceColor',[238,206,179]./255,'EdgeColor','none','FaceAlpha',0.6,'Parent',ha);
patch(struct('faces',mesh.tess_ctx.face,'vertices',mesh.tess_ctx.vert),'FaceColor',[1 0.7 0.7],'EdgeColor','none','Parent',ha);

view([az,el]);
axis(ha,'image','off')
material shiny
lighting gouraud
hl = camlight('headlight');
set(hf,'WindowButtonMotionFcn',@(~,~) camlight(hl,'headlight'),'Interruptible','off','busyaction','cancel');

% Sensors (coreg.)
h_sens  = plot3(D.sensors('MEG').chanpos(:,1),D.sensors('MEG').chanpos(:,2),D.sensors('MEG').chanpos(:,3),'og');
set(h_sens,'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize', 12,'MarkerEdgeColor','k');

for f = 1:3
    % MEG coils (coreg.)
    plot3(fid_polhemus(f,1),fid_polhemus(f,2),fid_polhemus(f,3),'+','MarkerSize',12,'color',cols{f});
    
    % MRI fiducials (coreg.)
   plot3(fid_mri(f,1),fid_mri(f,2),fid_mri(f,3),'o','MarkerSize',12,'color',cols{f});
end

% % MRI org fids
% plot3(fids.nas(1),fids.nas(2),fids.nas(3),'*','markersize',14,'Color',cols{1});
% plot3(fids.lpa(1),fids.lpa(2),fids.lpa(3),'*','markersize',14,'Color',cols{2});
% plot3(fids.rpa(1),fids.rpa(2),fids.rpa(3),'*','markersize',14,'Color',cols{3});

material dull % Improve visibility of the brain by reducing reflection glare

%indiv_mri_fid = [9.9898,142.5147,-8.1787;-55.7659,49.4636,-26.2089;88.7153,62.4787,-29.1394];
%spm_eeg_inv_transform_points



