clc; clear all; close all;

load leadfield
load elec_aligned
load tmp-headmodel
%%

leadfield_cell = leadfield.leadfield;
voxels_pos = leadfield.pos;
dim = leadfield.dim;
Ne = length(elec_aligned.label);
Nv = length(leadfield.pos);


%%
inside_sources = find(leadfield.inside);

real_indices = [500,900];

inside_index = inside_sources(real_indices);

% real_J_pos = voxels_pos(real_indices,:);
real_J = [1,0,0; -1,0,0];
voxel_J = zeros(Ne,0);
M = zeros(Ne,1);

for i= 1:length(real_indices)
    M = M + leadfield_cell{real_indices(i)}*real_J(i,:)';
end
%%
plot_haedmodel_elec_topo(headmodel.bnd(3),elec_aligned, M)
colorbar

for i = 1:length(inside_index)
    pos = [leadfield.pos(inside_index(i), 1);
        leadfield.pos(inside_index(i), 2);
        leadfield.pos(inside_index(i), 3)];
    hold on
    plot_real_dipole_pos(pos);
    plot_dipole_mom(pos, real_J(i,:))
    plot_dipole_number(pos,i)
    view(35,8)
end

%% 
index_cubic = reshape(1:prod(dim),dim);


%% Inverse problem

index = 500

x = 1 + mod(index,dim(1))
y = 



