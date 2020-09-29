load leadfield
%%

leadfield_cell = leadfield.leadfield; 
voxels_pos = leadfield.pos;
dim = leadfield.dim;
Ne = length(elec_aligned.label)
Nv = length(leadfield.pos);


%%
real_indices = [500,900];
real_J_pos = voxels_pos(real_indices,:);
real_J = [1,0,0; 1,0,0];
voxel_J = zeros(Ne,0);
M = zeros(Ne,1);

for i= 1:length(real_indices)
    M = M + leadfield_cell{real_indices(i)}*real_J(i,:)';
end
%%
hold on; ft_plot_topo3d(elec_aligned.elecpos, M)
hold on; ft_plot_sens(elec_aligned,'style','blue');
%     hold on; ft_plot_mesh(vol.bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on; ft_plot_mesh(headmodel.bnd(2),'edgecolor','none','facealpha',0.4);
hold on; ft_plot_mesh(headmodel.bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4]);

%%













