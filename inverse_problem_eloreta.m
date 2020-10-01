clc; clear all; close all;
addpath('utils');
load leadfield
load elec_aligned
load tmp-headmodel
%%

leadfield_cell = leadfield.leadfield;
voxels_pos = leadfield.pos;
dim = leadfield.dim;

inside_sources = find(leadfield.inside);
leadfield_inside = leadfield_cell(inside_sources);

Ne = length(elec_aligned.label);
Nv = length(leadfield_inside);

resolution = max(abs(diff(leadfield.pos(1:2,:))));

%%
real_indices = [500,900];

dipol_indices = inside_sources(real_indices);

% real_J_pos = voxels_pos(real_indices,:);
real_J_mom = [1,0,0; -1,0,0];
voxel_J = zeros(Ne,0);
M = zeros(Ne,1);

for i= 1:length(real_indices)
    M = M + leadfield_cell{real_indices(i)}*real_J_mom(i,:)';
end
%%
% plot_haedmodel_elec_topo(headmodel.bnd(3),elec_aligned, M)
% colorbar
hold on; ft_plot_mesh(headmodel.bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4],'facealpha',0.05);

plot_brain_grid_insides(leadfield,8)

for i = 1:length(dipol_indices)
    pos = [leadfield.pos(dipol_indices(i), 1);
        leadfield.pos(dipol_indices(i), 2);
        leadfield.pos(dipol_indices(i), 3)];
    hold on
    plot_real_dipole_pos(pos);
    plot_dipole_mom(pos, real_J_mom(i,:),resolution / 2)
    plot_dipole_number(pos,i)
    view(35,8)
end

%% 
index_cubic = reshape(1:prod(dim),dim);

%% Inverse problem -> loreta
% plot_brain_grid(leadfield)
clc
orien = 1;

K = zeros(Ne,Nv,1);
for i = 1:Ne
    K(:,i) = (leadfield_inside{i}(:,orien));
end


index_real = 500;
% leadfield_mat;
%%
W = eye(Nv,Nv);
H = eye(Ne) - ones(Ne,Ne)/Ne;
alpha = 0.01;
maxEpoch = 1000;
epsStop = 1e-5;
for epoch = 1:maxEpoch
    %         wi = [KT i (KW-1 KT + aH)+Ki]1/2
    W_save = W;
    C = pinv(K* pinv(W)* K' + alpha * H );
    for i = 1:Nv
        Ki = K(:,i);
        W(i,i) = (Ki' * C * Ki).^(0.5);
    end

    norm(W-W_save)
    if norm(W-W_save) < epsStop
        break
    end
end
J = zeros(Nv,1);
for i = 1:Nv
    Ki = K(:,i);
    J(i) = Ki'*((W(i,i)*(K*(W\K')+alpha*H))\M);
end

[~,index] = max(J)

err = norm(inside_pos(index,:) - inside_pos(index_real,:)) * 1000 %mm


%%
% c1 = 4; c2 = 4;
% for c = 1:length(dim(3))
%     subplot(c1,c2);
%     imagesc(
%     
%     
%     
%     
% end
    
%%% MohammadRaziei



