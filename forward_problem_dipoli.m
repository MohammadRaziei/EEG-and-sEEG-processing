clear; clc; close all;
addpath('utils')
%% load mri
mri = ft_read_mri('Subject01.mri');

%% segmentation

cfg=[];
cfg.output={'brain','skull','scalp'};
segmentedmri = ft_volumesegment(cfg,mri);

save tmp-segmentedmri segmentedmri;

%% triangulation
load tmp-segmentedmri;

cfg=[];
cfg.tissue={'brain','skull','scalp'};
% cfg.tissue    = {'gray','white','csf','skull','scalp'};
cfg.numvertices=[3000 2000 1000];
% cfg.numvertices=[2000 2000 2000 3000 2000 1000];
%cfg.sourceunits=segmentedmri.unit;
headmesh=ft_prepare_mesh(cfg,segmentedmri);
save tmp-headmesh headmesh;

clear segmentedmri;
%% headmodel
if(false)
% %% dipoli
load bnd3
cfg=[];
cfg.method='dipoli';
% cfg.conductivity = [0.33 0.14 1.79 0.01 0.43];           % order follows mesh.tissuelabel
cfg.conductivity = [0.75 0.01 0.43];           % order follows mesh.tissuelabel
headmodel=ft_prepare_headmodel(cfg,bnd);
save vol vol;

%%% bemcp
% load bnd3
% cfg=[];
% cfg.method='bemcp';
% cfg.conductivity = [0.75 0.01 0.43];
% vol =ft_prepare_headmodel(cfg,bnd);
% save vol3 vol

% %% Openmeeg
% % system('module load openmeeg'); % Load openmeeg paths
% cfg=[];
% cfg.method='openmeeg';
% vol_openmeeg=ft_prepare_headmodel(cfg,bnd);
% % warning- takes ~40 minutes
% save vol_openmeeg vol_openmeeg

end
load headmodel_dipoli_subject01
headmodel = vol; clear vol
save tmp-headmodel headmodel
%% visualization
% 
% figure;
% ft_plot_mesh(headmodel.bnd(1),'facecolor','none');  %skin
% figure;
% ft_plot_mesh(headmodel.bnd(2),'facecolor','none');
% figure;
% ft_plot_mesh(headmodel.bnd(3),'facecolor','none');

%%
figure;
ft_plot_mesh(headmodel.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(headmodel.bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(headmodel.bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4]);
view(50,0)
save_figure(gcf, "readme/forward-headmodel-mesh.png")

%% align electrodes

elec = ft_read_sens('standard_1020.elc');
%
% figure;
% ft_plot_sens(elec,'style','blue');
% hold on;
% ft_plot_mesh(vol.bnd(3),'facealpha', 0.85, 'edgecolor', 'none', 'facecolor', [0.65 0.65 0.65]); %scalp
% ;
% %
 nas=mri.hdr.fiducial.mri.nas;
 lpa=mri.hdr.fiducial.mri.lpa;
 rpa=mri.hdr.fiducial.mri.rpa;
%
 transm=mri.transform;
%
 nas=ft_warp_apply(transm,nas, 'homogenous');
 lpa=ft_warp_apply(transm,lpa, 'homogenous');
 rpa=ft_warp_apply(transm,rpa, 'homogenous');
%
%  fiducials.chanpos=[nas; lpa; rpa];
%  fiducials.label={'Nz','LPA','RPA'};
%  fiducials.unit='mm';
%
%  % ensure the elec to have a coordsys field, in order to avoid
%  % the interactive step due to the call to ft_determine_coordsys
%  % in ft_electroderealign
%  elec.coordsys = 'ctf';
%
%  cfg=[];
%  cfg.method='fiducial';
%  cfg.template=fiducials;
% %  cfg.elec = elec;
%  cfg.fiducial={'Nz', 'LPA', 'RPA'};
%  elec_align=ft_electroderealign(cfg, elec);


% create a structure similar to a template set of electrodes
fid.elecpos       = [nas; lpa; rpa];       % ctf-coordinates of fiducials
fid.label         = {'Nz','LPA','RPA'};    % same labels as in elec
fid.unit          = 'mm';                  % same units as mri

% alignment
cfg               = [];
cfg.method        = 'fiducial';
cfg.target        = fid;                   % see above
cfg.elec          = elec;
cfg.fiducial      = {'Nz', 'LPA', 'RPA'};  % labels of fiducials in fid and in elec

elec_align        = ft_electroderealign(cfg);

save elec_align elec_align;


%%
% %   cfg=[];
% %   cfg.method='interactive';
% %   cfg.elec=elec_align;
% %   cfg.headshape=vol.bnd(1);
% % %
% %   elec_align=ft_electroderealign(cfg);
% %  ;
%
%
%
%   figure;
% ft_plot_sens(elec_align,'style','sk');
% hold on;
% ft_plot_mesh(vol.bnd(1),'facealpha', 0.85, 'edgecolor', 'none', 'facecolor', [0.65 0.65 0.65]); %scalp
%
%   ;
% %% concentric spheres
% load bnd;
% cfg=[];
% cfg.method='concentricspheres';
% cfg.conductivity=[0.33 0.041 0.33];
% vol_cs=ft_prepare_headmodel(cfg,bnd);
%
% figure;
% ft_plot_vol(vol_cs,'facealpha', 0.3)
% hold on;
% ft_plot_mesh(bnd(3),'facealpha', 0.3, 'facecolor', 'red', 'edgecolor', 'none');
%
% figure;
% ft_plot_vol(vol_cs,'facealpha', 0.3)
% hold on;
% ft_plot_mesh(bnd(1),'facealpha', 0.3, 'facecolor', 'red', 'edgecolor', 'none');
%
%%
% calculate BEM leadfield
% for ll=1:length(vertic)
%   cfg = [];
%   cfg.sourcemodel = sourcemodel;
%   cfg.headmodel = volbem{ll};
%   cfg.elec = elec;
%   leadfieldBEM{ll} = ft_prepare_leadfield(cfg);
% end




transform = [
    1 0 0 15
    0 1 0 0
    0 0 1 -10
    0 0 0 1
    ];
elec_aligned = ft_transform_sens(transform, elec_align);

save elec_aligned elec_aligned
%%%%%%%%%
figure;
ft_plot_sens(elec_aligned,'style','blue');
hold on;
ft_plot_mesh(headmodel.bnd(3),'facealpha', 0.85, 'edgecolor', 'none', 'facecolor', [0.65 0.65 0.65]); %scalp
hold on;
ft_plot_mesh(headmodel.bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(headmodel.bnd(1),'edgecolor','none','facecolor',[0.4 0.6 0.4]);
%%
close all;
%

cfg = [];
cfg.vol = headmodel;
cfg.grid.resolution = 10;
cfg.elec = elec_aligned;
cfg.grid.unit = 'mm';
cfg.inwardshift = 5;
sourcemodel = ft_prepare_sourcemodel(cfg);
save sourcemodel sourcemodel
%% calculate theoretical single sphere leadfield
cfg = [];
cfg.grid = sourcemodel;
cfg.headmodel = headmodel;
cfg.elec = elec_aligned;
leadfield = ft_prepare_leadfield(cfg);
save leadfield leadfield



%%
% index = 1000;
% f = 5;
% figure;
% for t = 0:0.01:0
%     q_mom = [1 0 0]';
%     q_pos = leadfield.pos(index);
%     q = q_mom * cos(2*pi*f*t);
%     
%     G = leadfield.leadfield{index};
%     M = G*q;
%     
%     clf
%     hold on; ft_plot_topo3d(elec.elecpos, M)
%     hold on; ft_plot_sens(elec_aligned,'style','blue');
% %     hold on; ft_plot_mesh(vol.bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
%     hold on; ft_plot_mesh(headmodel.bnd(2),'edgecolor','none','facealpha',0.4);
%     hold on; ft_plot_mesh(headmodel.bnd(1),'edgecolor','none','facecolor',[0.4 0.6 0.4]);
%     
%     pause(0.001)
% end